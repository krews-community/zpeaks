package step.subpeaks

import model.*
import util.*
import java.lang.Math.*

/**
 * Represents a Gaussian distribution with shape (adds skewness).
 * Comparison functions only compare amplitude (used for sorting
 * from most to least important when curve fitting).
 *
 * shape: shape parameter (alpha) capturing skewness.
 * amplitude: suqare root of the amplitude of the distribution.
 * mean: center point of the distribution.
 * stdev: standard deviation of the distribution.
 */
data class SkewGaussianParameters(val amplitude: Double, val mean: Double, val stdDev: Double, val shape: Double)

data class StandardGaussianParameters(val amplitude: Double, val mean: Double, val stdDev: Double)

interface GaussianParameters

data class PeakGaussian<P : GaussianParameters>(val region: Region, val gaussianParameters: P)

/**
 * Performs a gaussian step.subpeaks.fit on the given region and pushes output to the passed vectors
 * for scores, regions, and parameters. This method will retry the step.subpeaks.fit up to five times
 * for large errors (>0.1) and splits large regions (>5kb) at the minimum element to avoid
 * very long step.subpeaks.fit times.
 *
 * values: the vector of values to step.subpeaks.fit.
 * region: the coordinates of the region being step.subpeaks.fit.
 * tregionout: output vector receiving coordinates of step.subpeaks.fit subpeaks.
 * tscoreout: output vector receiving region scores.
 * gaussiansout: output vector receiving parameters of step.subpeaks.fit subpeaks.
 */
fun <P> fit(pileUpValues: List<Int>, pileUpRegion: Region): List<P> {
    // Subtract out background
    val min = pileUpValues.min()!!
    val valuesWithoutBackground = pileUpValues.map { it - min }

    // Put all on same scale
    val avg = valuesWithoutBackground.average()
    val scaledValues = valuesWithoutBackground.map { it / avg }

    // Get initial candidates list
    val candidates = scaleSpaceZeros(scaledValues)

    /*
     * here we compute guesses for the Gaussians' parameters based on zero crossing locations
     * the means and stddevs are computed with simple equations
     * the amplitudes are fit to the data using least squares after the means and stddevs are computed
     */
    val gaussians =

}

/*
 * Uses the scale-space image of the data to generate an initial list
 * of candidate Gaussians. Zero crossings as they emerge from most blurred
 * to least define the Gaussians' parameters.
 *
 * values: input vector from which to generate the scale space image.
 */
private fun scaleSpaceZeros(values: List<Double>): List<Region> {
    val candidates = mutableListOf<Region>()

    // Start with a kernel near the square root of the vector size and loop down to 2
    val kernel = Math.floor(sqrt(values.size.toDouble()) * 1.1).toInt()
    for (kernelSize in kernel downTo 2 step 1) {

        // Smooth at this kernel width
        val blurred = scaleSpaceSmooth(values, kernelSize.toDouble())

        // Find zero crossings
        val zeroCrossings = mutableListOf<Int>()
        var cSign = blurred[0] > 0
        for (j in 1 until blurred.size) {
            if (cSign != blurred[j] > 0) {
                cSign = blurred[j] > 0
                zeroCrossings += j
            }
        }

        // Match existing to nearest zero crossings in current iteration
        val matched = MutableList(zeroCrossings.size) { false }
        val candidateIter = candidates.listIterator()
        candidateIter.forEach { candidate ->
            // Find match for start and stop
            var cStartMin = 1E308.toInt()
            var cStopMin = 1E308.toInt()
            var startMatch = -1
            var stopMatch = -1
            for (crossingIndex in 0 until zeroCrossings.size) {
                if (matched[crossingIndex]) continue
                val startDist = abs(candidate.start - crossingIndex)
                val stopDist = abs(candidate.end - crossingIndex)
                val len = candidate.end - candidate.start
                if (startDist < cStartMin && startDist < len) {
                    startMatch = crossingIndex
                    cStartMin = startDist
                } else if (stopDist <= cStopMin && stopDist < len) {
                    stopMatch = crossingIndex
                    cStopMin = stopDist
                }
            }

            // Skip if no match else adjust
            if (stopMatch != -1 && startMatch != -1) {
                matched[startMatch] = true
                matched[stopMatch] = true
                candidateIter.set(Region(zeroCrossings[startMatch], zeroCrossings[stopMatch]))
            }
        }

        // append unmatched zero crossings as new candidates
        val unmatched = mutableListOf<Int>()
        for ((crossingIndex, crossing) in zeroCrossings.withIndex()) {
            if (!matched[crossingIndex]) unmatched += crossing
        }
        for (n in 1 until unmatched.size step 2) {
            candidates += Region(unmatched[n - 1], unmatched[n])
        }
    }

    // Remove any uncompleted candidates before return
    candidates.removeAll { c -> c.start == -1 || c.end == -1 }

    return candidates
}

fun splitAtMin(values: List<Double>): List<List<Double>> {
    // look only within the middle half to make sure we shrink the sizes by a reasonable amount
    val quarterSize = values.size / 4
    val minIndex = values.subList(quarterSize, quarterSize * 3).indexOfMin()!! + quarterSize

    // Return two halves, split at min index.
    return listOf(values.subList(0, minIndex), values.subList(minIndex, values.size))
}