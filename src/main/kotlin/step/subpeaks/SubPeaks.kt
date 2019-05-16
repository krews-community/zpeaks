package step.subpeaks

import model.*
import mu.KotlinLogging
import org.apache.commons.math3.linear.*
import org.apache.commons.math3.linear.MatrixUtils.*
import org.apache.commons.math3.util.FastMath.*
import util.*
import java.lang.Double.isNaN
import java.lang.Exception

private val log = KotlinLogging.logger {}

interface GaussianParameters {
    var amplitude: Double
    var mean: Double
    var stdDev: Double
}

data class CandidateGaussian<T : GaussianParameters>(
    val region: Region,
    val parameters: T
)

data class Fit<T: GaussianParameters>(
    val subPeaks: List<SubPeak<T>>,
    val error: Double
)

data class SubPeak<T : GaussianParameters>(
    val region: Region,
    val score: Double,
    val gaussianParameters: T
)

data class OptimizeResults<T : GaussianParameters>(
    val parameters: List<T>,
    val error: Double,
    val iterations: Int
)

typealias Optimizer<T> = (values: DoubleArray, candidateGaussians: List<CandidateGaussian<T>>, gaussians: List<T>) -> OptimizeResults<T>

const val SOFT_MAX_CANDIDATES = 5
const val SOFT_MAX_CUT_RATIO = 0.05
const val HARD_MAX_CANDIDATES = 15

class SubPeakFitter<T: GaussianParameters> (
    private val initParameters: (region: Region) -> T,
    private val optimize: Optimizer<T>,
    private val parametersToRegion: (parameters: T, offset: Int) -> Region
) {

    /**
     * Performs a gaussian fit on the given region and pushes output to the passed vectors
     * for scores, regions, and parameters.
     */
    fun fit(values: List<Double>, start: Int): Fit<T> {
        // Subtract out background
        val min = values.min()!!
        val valuesWithoutBackground = values.map { it - min }

        // Get initial candidates list
        val candidates = findCandidates(valuesWithoutBackground)

        val splitPeaks = split(PreparedPeakData(start, valuesWithoutBackground, candidates))
        val fits = splitPeaks.map { fitCandidates(it) }

        return Fit(
            subPeaks = fits.flatMap { it.subPeaks },
            error = fits.sumByDouble { it.error } / fits.size
        )
    }

    private data class PreparedPeakData(
        val start: Int,
        val valuesWithoutBackground: List<Double>,
        val candidates: List<Region>
    )

    private fun split(peak: PreparedPeakData): List<PreparedPeakData> {
        val splitPeaks = mutableListOf<PreparedPeakData>()
        var peaksToSplit = mutableListOf(peak)

        while (peaksToSplit.isNotEmpty()) {
            val nextPeaksToSplit = mutableListOf<PreparedPeakData>()
            for (p in peaksToSplit) {
                if (shouldSplit(peak)) {
                    val minIndex = p.valuesWithoutBackground.indexOfMin()!!
                    nextPeaksToSplit += PreparedPeakData(
                        start = p.start,
                        valuesWithoutBackground = p.valuesWithoutBackground.subList(0, minIndex),
                        candidates = p.candidates.filter { it.start < minIndex }
                    )
                    nextPeaksToSplit += PreparedPeakData(
                        start = p.start + minIndex,
                        valuesWithoutBackground = p.valuesWithoutBackground.subList(minIndex, p.valuesWithoutBackground.size),
                        candidates = p.candidates.filter { it.start > minIndex }
                    )
                } else {
                    splitPeaks += p
                }
            }
            peaksToSplit = nextPeaksToSplit
        }

        return splitPeaks.sortedBy { it.start }
    }

    private fun shouldSplit(peak: PreparedPeakData): Boolean {
        if (peak.candidates.size <= SOFT_MAX_CANDIDATES) return false
        if (peak.candidates.size > HARD_MAX_CANDIDATES) return true
        val minMaxRatio = peak.valuesWithoutBackground.min()!! / peak.valuesWithoutBackground.max()!!
        return minMaxRatio <= SOFT_MAX_CUT_RATIO
    }

    private fun bestSplitIndex(values: List<Double>, candidates: List<Region>) {
        //val localMinima = values.localMinima()
        TODO("")
    }

    private fun fitCandidates(peak: PreparedPeakData): Fit<T> {
        // Put all on same scale
        val avg = peak.valuesWithoutBackground.average()
        val scaledValues = peak.valuesWithoutBackground.map { it / avg }.toDoubleArray()

        /*
         * here we compute guesses for the Gaussians' parameters based on zero crossing locations
         * the means and standard deviations are computed with simple equations
         * the amplitudes are fit to the data using least squares after the means and standard deviations are computed
         */
        val candidateGaussians = peak.candidates.map { CandidateGaussian(it, initParameters(it)) } as MutableList

        /*
         * Solve the system of linear equations to get the amplitude guesses
         */
        val lim = peak.valuesWithoutBackground.size
        val m = Array(candidateGaussians.size) { DoubleArray(candidateGaussians.size) }
        val d = DoubleArray(candidateGaussians.size)
        val distributions = candidateGaussians
            .map { gaussianDistribution(1.0, it.parameters.mean, it.parameters.stdDev, lim) }

        // Set coefficients: sum all bp x of (sum all gaussians j of (kth gaussian(x) * jth gaussian(x)))
        for (j in 0 until candidateGaussians.size) {
            for (k in 0 until candidateGaussians.size) {
                var value = 0.0
                for (n in 0 until lim) {
                    value += distributions[k][n] * distributions[j][n]
                }
                m[k][j] = value
            }
        }

        // Set values: sum all bp x of (curve value x * jth gaussian(x))
        for (j in 0 until candidateGaussians.size) {
            var value = 0.0
            for (n in 0 until lim) {
                value += distributions[j][n] * scaledValues[n]
            }
            d[j] = value
        }

        val r = LUDecomposition(createRealMatrix(m)).solver.solve(createRealVector(d))

        for (j in 0 until candidateGaussians.size) {
            val s = if (r.getEntry(j) < 0) -1.0 else 1.0
            // Set Amplitude
            candidateGaussians[j].parameters.amplitude = sqrt(abs(r.getEntry(j)) * candidateGaussians[j].parameters.stdDev) * s
        }
        // Sort gaussians by amplitude
        candidateGaussians.sortByDescending { it.parameters.amplitude * it.parameters.amplitude / it.parameters.stdDev }

        log.debug { "Candidate gaussians size: ${candidateGaussians.size}" }
        log.debug { "Candidate gaussians: $candidateGaussians" }

        /*
         * Perform the actual fit of the curve to a sum of gaussians.size Gaussians
         * optimization is performed with the Levenberg-Marquardt algorithm
         */
        var bestOptimized: OptimizeResults<T>? = null
        var previousOptimized: OptimizeResults<T>? = null
        for (j in 0 until candidateGaussians.size) {
            val currentGuess =
                candidateGaussians.subList(0, j+1).map { it.parameters }
            //(previousOptimized?.parameters ?: listOf()) + candidateGaussians[j].parameters

            val optimizeResults = optimize(scaledValues, candidateGaussians.subList(0, j+1), currentGuess)
            log.debug { "Attempt #$j. Current guess size = ${currentGuess.size}. Error = ${optimizeResults.error}" }

            previousOptimized = optimizeResults
            if (optimizeResults.error < 0.05) {
                log.debug { "Answer has acceptable error < 0.05!" }
                bestOptimized = optimizeResults
                break
            }
            if (bestOptimized == null || optimizeResults.error < bestOptimized.error) {
                bestOptimized = optimizeResults
            }
        }

        log.debug { "Optimized Result size: ${bestOptimized!!.parameters.size}" }
        log.debug { "Optimized Result error: ${bestOptimized!!.error}" }
        log.debug { "Optimized Result iterations: ${bestOptimized!!.iterations}" }
        log.debug { "Optimized Result: $bestOptimized" }
        val sqrtAvg = sqrt(avg)
        // Bring parameters back to original scale and add background back in. Also add offset to mean
        bestOptimized!!.parameters.forEach {
            it.amplitude *= sqrtAvg
            it.mean += peak.start
        }

        val subPeaks = bestOptimized.parameters
            .sortedBy { it.mean }
            .map { SubPeak(parametersToRegion(it, peak.start), parametersToScore(it), it) }

        return Fit(subPeaks, bestOptimized.error)
    }
}

fun parametersToScore(parameters: GaussianParameters): Double = parameters.amplitude.pow(2) / parameters.stdDev

/*
 * Uses the scale-space image of the data to generate an initial list
 * of candidate Gaussians. Zero crossings as they emerge from most blurred
 * to least define the Gaussians' parameters.
 *
 * values: input vector from which to generate the scale space image.
 */
fun findCandidates(values: List<Double>): List<Region> {
    // This method is slow for array with length > 5000
    if (values.size > 5000) {
        val (first, second) = splitAtMin(values)
        return findCandidates(first) + findCandidates(second)
    }

    var candidates = mutableListOf<Region>()

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

        // Find matching zero-crossings in existing candidates
        val currentCandidates = mutableListOf<Region>()
        val matchedZeroCrossings = mutableSetOf<Int>()
        for (candidate in candidates) {
            val candidateLength = candidate.end - candidate.start
            // Find the closest matching zero crossings within candidateLength to the current candidate
            var startMatch: Int? = null
            var startMatchDist: Int? = null
            var endMatch: Int? = null
            var endMatchDist: Int? = null
            for (crossing in zeroCrossings) {
                if (matchedZeroCrossings.contains(crossing)) continue
                val startDist = abs(crossing - candidate.start)
                val endDist = abs(crossing - candidate.end)
                if (startDist < candidateLength && (startMatch == null || startDist < startMatchDist!!)) {
                    startMatch = crossing
                    startMatchDist = startDist
                } else if (endDist < candidateLength && (endMatch == null || endDist < endMatchDist!!)) {
                    endMatch = crossing
                    endMatchDist = endDist
                }
            }

            // If we got a match, add new match to candidates in place of old. Otherwise add old candidate.
            if (startMatch != null && endMatch != null) {
                matchedZeroCrossings += startMatch
                matchedZeroCrossings += endMatch
                currentCandidates += Region(startMatch, endMatch)
            } else {
                currentCandidates += candidate
            }
        }

        // Add unmatched zero crossings as new candidates
        val unmatchedZeroCrossings = zeroCrossings.filter { !matchedZeroCrossings.contains(it) }
        for (n in 1 until unmatchedZeroCrossings.size step 2) {
            currentCandidates += Region(unmatchedZeroCrossings[n - 1], unmatchedZeroCrossings[n])
        }

        candidates = currentCandidates
    }

    return candidates
}
