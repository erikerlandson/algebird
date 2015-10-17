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

import org.scalatest._

import com.twitter.algebird.matchers.seq._

class TDigestTest extends FlatSpec with Matchers {
  import org.apache.commons.math3.distribution.RealDistribution

  val seed = 235711L
  scala.util.Random.setSeed(seed)

  val ss = 100000
  val delta = 50.0 / 1000

  val maxD = 0.05
  val maxDI = 0.05

  def testDistribution(dist: RealDistribution, stdv: Double): Boolean = {
    import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest

    dist.reseedRandomGenerator(seed)

    val td = TDigest.sketch(Iterator.fill(ss) { dist.sample }, delta = delta)

    val xmin = td.clusters.keyMin.get
    val xmax = td.clusters.keyMax.get
    val step = (xmax - xmin) / 1000
    val d = (xmin to xmax by step).iterator
      .map(x => math.abs(td.cdf(x) - dist.cumulativeProbability(x))).max

    val dInv = (0.01 to 0.99 by 0.01).iterator
      .map(x => math.abs(td.cdfInverse(x) - dist.inverseCumulativeProbability(x))).max / stdv

    val pass = d <= maxD && dInv <= maxDI
    if (!pass) println(s"d= $d  dInv= $dInv")
    pass
  }

  it should "sketch a uniform distribution" in {
    import org.apache.commons.math3.distribution.UniformRealDistribution
    val dist = new UniformRealDistribution()
    testDistribution(dist, math.sqrt(dist.getNumericalVariance())) should be (true)
  }

  it should "sketch a normal distribution" in {
    import org.apache.commons.math3.distribution.NormalDistribution
    val dist = new NormalDistribution()
    testDistribution(dist, math.sqrt(dist.getNumericalVariance())) should be (true)
  }

  it should "sketch an exponential distribution" in {
    import org.apache.commons.math3.distribution.ExponentialDistribution
    val dist = new ExponentialDistribution(1.0)
    testDistribution(dist, math.sqrt(dist.getNumericalVariance())) should be (true)
  }
}
