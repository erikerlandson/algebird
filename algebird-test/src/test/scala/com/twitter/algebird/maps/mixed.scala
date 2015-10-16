/*
Copyright 2012 Twitter, Inc.

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

package com.twitter.algebird.maps

import org.scalatest._

import com.twitter.algebird.Monoid

import com.twitter.algebird.matchers.seq._

object mixed {
  import math.Numeric
  import com.twitter.algebird.maps.increment._
  import com.twitter.algebird.maps.prefixsum._
  import com.twitter.algebird.maps.nearest._

  object tree {
    import com.twitter.algebird.maps.increment.tree._
    import com.twitter.algebird.maps.prefixsum.tree._
    import com.twitter.algebird.maps.nearest.tree._

    trait NodeMix[K, V, P] extends NodePS[K, V, P] with NodeInc[K, V] with NodeNearMap[K, V]

    trait LNodeMix[K, V, P] extends NodeMix[K, V, P]
      with LNodePS[K, V, P] with LNodeInc[K, V] with LNodeNearMap[K, V]

    trait INodeMix[K, V, P] extends NodeMix[K, V, P]
      with INodePS[K, V, P] with INodeInc[K, V] with INodeNearMap[K, V] {
      val lsub: NodeMix[K, V, P]
      val rsub: NodeMix[K, V, P]
    }
  }

  import tree._

  object infra {
    import com.twitter.algebird.maps.redblack.tree._
    import com.twitter.algebird.maps.ordered.tree.DataMap

    class Inject[K, V, P](
      val keyOrdering: Numeric[K],
      val valueMonoid: Monoid[V],
      val prefixMonoid: IncrementingMonoid[P, V]) {

      def iNode(clr: Color, dat: Data[K], ls: Node[K], rs: Node[K]) =
        new Inject[K, V, P](keyOrdering, valueMonoid, prefixMonoid) with INodeMix[K, V, P] with MixedMap[K, V, P] {
          // INode[K]
          val color = clr
          val lsub = ls.asInstanceOf[NodeMix[K, V, P]]
          val rsub = rs.asInstanceOf[NodeMix[K, V, P]]
          val data = dat.asInstanceOf[DataMap[K, V]]
          // INodePS[K, V, P]
          val prefix = prefixMonoid.inc(prefixMonoid.plus(lsub.pfs, rsub.pfs), data.value)
          // INodeNear[K, V]
          val kmin = lsub match {
            case n: INodeMix[K, V, P] => n.kmin
            case _ => data.key
          }
          val kmax = rsub match {
            case n: INodeMix[K, V, P] => n.kmax
            case _ => data.key
          }
        }
    }

  }

  import infra._

  sealed trait MixedMap[K, V, P]
    extends IncrementMapLike[K, V, INodeMix[K, V, P], MixedMap[K, V, P]]
    with PrefixSumMapLike[K, V, P, INodeMix[K, V, P], MixedMap[K, V, P]]
    with NearestMapLike[K, V, INodeMix[K, V, P], MixedMap[K, V, P]] {

    override def toString =
      "MixedMap(" +
        iterator.zip(prefixSumsIterator())
        .map(x => s"${x._1._1} -> (${x._1._2}, ${x._2})").mkString(", ") +
        ")"
  }

  object MixedMap {
    def key[K](implicit num: Numeric[K]) = new AnyRef {
      def value[V](implicit vm: Monoid[V]) = new AnyRef {
        def prefix[P](implicit im: IncrementingMonoid[P, V]): MixedMap[K, V, P] =
          new Inject[K, V, P](num, vm, im) with LNodeMix[K, V, P] with MixedMap[K, V, P]
      }
    }
  }
}

class MixedMapSpec extends FlatSpec with Matchers {
  import scala.language.reflectiveCalls

  import com.twitter.algebird.maps.prefixsum.IncrementingMonoid

  import com.twitter.algebird.maps.ordered.RBProperties._
  import com.twitter.algebird.maps.ordered.OrderedMapProperties._
  import com.twitter.algebird.maps.prefixsum.PrefixSumMapProperties._
  import com.twitter.algebird.maps.increment.IncrementMapProperties._
  import com.twitter.algebird.maps.nearest.NearestMapProperties._

  import mixed.MixedMap

  def mapType1 =
    MixedMap.key[Double].value[Int]
      .prefix(IncrementingMonoid.fromMonoid(implicitly[Monoid[Int]]))

  it should "pass randomized tree patterns" in {
    val data = Vector.tabulate(50)(j => (j.toDouble, j))
    (1 to 1000).foreach { u =>
      val shuffled = scala.util.Random.shuffle(data)
      val map = shuffled.foldLeft(mapType1)((m, e) => m + e)

      testRB(map)
      testKV(data, map)
      testDel(data, map)
      testPrefix(data, map)
      testIncrement(data, map)
      testNearest(data, map)
    }
  }
}
