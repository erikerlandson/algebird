/*
Copyright 2015 Twitter, Inc.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

package com.twitter.algebird

import tdmap.TDigestMap

/** A t-digest object */
case class TDigest(
  delta: Double,
  recluster: Int,
  nclusters: Int,
  clusters: TDigestMap) {

  private case class Cluster(centroid: Double, mass: Double, massUB: Double)

  /**
   * Returns a new t-digest with value x included in its sketch; td + x is equivalent to
   * td + (x, 1).
   */
  def +[N](x: N)(implicit num: Numeric[N]): TDigest = this.+((x, 1))

  /**
   * Returns a new t-digest with new pair (x, w) included in its sketch.
   * This implements 'algorithm 1' from:
   * Computing Extremely Accurate Quantiles Using t-Digests
   * Ted Dunning and Otmar Ertl
   * https://github.com/tdunning/t-digest/blob/master/docs/t-digest-paper/histo.pdf
   */
  def +[N1, N2](xw: (N1, N2))(implicit num1: Numeric[N1], num2: Numeric[N2]): TDigest = {
    val s = this.update(xw)
    if (s.nclusters < recluster) s
    else {
      // too many clusters: attempt to compress it by re-clustering
      val ds = scala.util.Random.shuffle(s.clusters.toVector)
      ds.foldLeft(TDigest.empty(delta, recluster))((d, e) => d.update(e))
    }
  }

  private def update[N1, N2](xw: (N1, N2))(implicit num1: Numeric[N1], num2: Numeric[N2]) = {
    val xn = num1.toDouble(xw._1)
    var wn = num2.toDouble(xw._2)
    require(wn > 0.0, "data weight must be > 0")

    // Get the current cluster nearest to incoming (xn)
    // Note: 'near' will have length 0,1, or 2:
    // length 0 => current cluster map was empty (no data yet)
    // length 1 => exactly one cluster was closest to (xn)
    // length 2 => (xn) was mid-point between two clusters (both are returned, in key order)
    val near = clusters.nearest(xn)

    if (near.isEmpty) {
      // our map is empty, so insert this pair as the first cluster
      TDigest(delta, recluster, nclusters + 1, clusters + ((xn, wn)))
    } else {
      // compute upper bounds for cluster masses, from their quantile estimates
      var massPS = clusters.prefixSum(near.head._1, open = true)
      val massTotal = clusters.sum
      val s = near.map {
        case (c, m) =>
          val q = (massPS + m / 2.0) / massTotal
          val ub = massTotal * delta * q * (1.0 - q)
          massPS += m
          Cluster(c, m, ub)
      }

      // assign new mass (wn) among the clusters
      var cmNew = Vector.empty[(Double, Double)]
      scala.util.Random.shuffle(s).foreach { clust =>
        if (wn <= 0.0) {
          // if we have already distributed all the mass, remaining clusters unchanged
          cmNew = cmNew :+ ((clust.centroid, clust.mass))
        } else if (xn == clust.centroid) {
          // if xn lies exactly on the centroid, add all mass in regardless of bound
          cmNew = cmNew :+ ((clust.centroid, clust.mass + wn))
          wn = 0.0
        } else if (clust.mass < clust.massUB) {
          // cluster can accept more mass, respecting its upper bound
          val dm = math.min(wn, clust.massUB - clust.mass)
          val mass = clust.mass + dm
          val dc = dm * (xn - clust.centroid) / mass
          wn -= dm
          cmNew = cmNew :+ ((clust.centroid + dc, mass))
        } else {
          // cluster is at its upper bound for mass, it remains unchanged
          cmNew = cmNew :+ ((clust.centroid, clust.mass))
        }
      }

      // any remaining mass becomes a new cluster
      if (wn > 0.0) cmNew = cmNew :+ ((xn, wn))

      // remove original clusters and replace with the new ones
      val clustDel = near.iterator.map(_._1).foldLeft(clusters)((c, e) => c - e)
      val clustNew = cmNew.foldLeft(clustDel)((c, p) => c.increment(p._1, p._2))
      val nc = nclusters - s.length + cmNew.length

      // return the updated t-digest
      TDigest(delta, recluster, nc, clustNew)
    }
  }

  def cdf[N](xx: N)(implicit num: Numeric[N]) = clusters.cdf(xx)
  def cdfInverse[N](qq: N)(implicit num: Numeric[N]) = clusters.cdfInverse(qq)
}

object TDigest {
  import scala.language.reflectiveCalls
  import com.twitter.algebird.Monoid

  /** return an empty t-digest */
  def empty(
    delta: Double = 0.5,
    recluster: Int = 1000) = {
    require(delta > 0.0, s"delta was not > 0")
    TDigest(delta, recluster, 0, TDigestMap.empty)
  }

  /** return a t-digest constructed from some data */
  def sketch[N](
    data: TraversableOnce[N],
    delta: Double = 0.5,
    recluster: Int = 1000)(implicit num: Numeric[N]) = {
    require(delta > 0.0, s"delta was not > 0")
    val td = data.foldLeft(empty(delta, recluster))((c, e) => c + ((e, 1)))
    val ds = scala.util.Random.shuffle(td.clusters.toVector)
    ds.foldLeft(empty(delta, recluster))((c, e) => c + e)
  }
}
