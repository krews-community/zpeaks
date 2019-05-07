package step.subpeaks

import model.*
import kotlin.math.*

data class StandardGaussianParameters(
    override var amplitude: Double,
    override var mean: Double,
    override var stdDev: Double
) : GaussianParameters

fun initStandardParameters(region: Region) = StandardGaussianParameters(
    amplitude = 0.0,
    mean = (region.start + region.end) / 2.0,
    stdDev = (region.end - region.start) / 2.0
)

fun curveValue(amplitude: Double, mean: Double, stdDev: Double, x: Double): Double =
        amplitude.pow(2) / stdDev * exp( -(x - mean) * (x - mean) / stdDev.pow(2) / 2)