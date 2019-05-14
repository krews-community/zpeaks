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

typealias Optimizer<T> = (values: DoubleArray, gaussians: List<T>, lambda: Double) -> OptimizeResults<T>

/**
 * Performs a gaussian fit on the given region and pushes output to the passed vectors
 * for scores, regions, and parameters. This method will retry the fit up to five times
 * for large errors (>0.1) and splits large regions (>5kb) at the minimum element to avoid
 * very long fit times.
 *
 * values: the vector of values to fit.
 * pileUpRegion: the coordinates of the region being fit.
 */
fun <T : GaussianParameters> fit(
    values: List<Double>,
    start: Int,
    initParameters: (region: Region) -> T,
    optimize: Optimizer<T>,
    parametersToRegion: (parameters: T, offset: Int) -> Region
): List<SubPeak<T>> {
    if (values.size > 5000) {
        val splitValues = splitAtMin(values)
        log.info { "Splitting: ${values.size} to ${splitValues.first.size} and ${splitValues.second.size}" }
        val firstFit = fit(splitValues.first, start, initParameters, optimize, parametersToRegion)
        val secondFit = fit(
            splitValues.second, start + splitValues.first.size,
            initParameters, optimize, parametersToRegion
        )
        return firstFit + secondFit
    }

    // Subtract out background
    val min = values.min()!!
    val valuesWithoutBackground = values.map { it - min }

    // Put all on same scale
    val avg = valuesWithoutBackground.average()
    val scaledValues = valuesWithoutBackground.map { it / avg }.toDoubleArray()

    // Get initial candidates list
    val candidates = findCandidates(scaledValues)

    /*
     * here we compute guesses for the Gaussians' parameters based on zero crossing locations
     * the means and standard deviations are computed with simple equations
     * the amplitudes are fit to the data using least squares after the means and standard deviations are computed
     */
    val gaussians = candidates.map { initParameters(it) } as MutableList
    log.debug { "Candidate gaussians size: ${gaussians.size}" }
    log.debug { "Candidate gaussians: $gaussians" }

    /*
     * Solve the system of linear equations to get the amplitude guesses
     */
    val lim = values.size
    val m = Array(gaussians.size) { DoubleArray(gaussians.size) }
    val d = DoubleArray(gaussians.size)
    val distributions = gaussians.map { gaussianDistribution(1.0, it.mean, it.stdDev, lim) }

    // Set coefficients: sum all bp x of (sum all gaussians j of (kth gaussian(x) * jth gaussian(x)))
    for (j in 0 until gaussians.size) {
        for (k in 0 until gaussians.size) {
            var value = 0.0
            for (n in 0 until lim) {
                value += distributions[k][n] * distributions[j][n]
            }
            m[k][j] = value
        }
    }

    // Set values: sum all bp x of (curve value x * jth gaussian(x))
    for (j in 0 until gaussians.size) {
        var value = 0.0
        for (n in 0 until lim) {
            value += distributions[j][n] * scaledValues[n]
        }
        d[j] = value
    }

    val r = LUDecomposition(createRealMatrix(m)).solver.solve(createRealVector(d))

    for (j in 0 until gaussians.size) {
        val s = if (r.getEntry(j) < 0) -1.0 else 1.0
        // Set Amplitude
        gaussians[j].amplitude = sqrt(abs(r.getEntry(j)) * gaussians[j].stdDev) * s
    }
    // Sort gaussians by amplitude
    gaussians.sortBy { it.amplitude * it.amplitude / it.stdDev }

    /*
     * Perform the actual fit of the curve to a sum of gaussians.size Gaussians
     * optimization is performed with the Levenberg-Marquardt algorithm
     */
    var bestOptimized: OptimizeResults<T>? = null
    var previousOptimized: OptimizeResults<T>? = null
    for (j in 0 until gaussians.size) {
        val currentGuess = if (previousOptimized != null) {
            previousOptimized.parameters + gaussians[j]
        } else {
            gaussians.subList(0, j+1)
        }

        val optimizeResults = optimize(scaledValues, currentGuess, 0.1)
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
        it.mean += start
    }

    return bestOptimized.parameters.map { SubPeak(parametersToRegion(it, start), parametersToScore(it), it) }
}

fun parametersToScore(parameters: GaussianParameters): Double = parameters.amplitude.pow(2) / parameters.stdDev

/*
 * Uses the scale-space image of the data to generate an initial list
 * of candidate Gaussians. Zero crossings as they emerge from most blurred
 * to least define the Gaussians' parameters.
 *
 * values: input vector from which to generate the scale space image.
 */
fun findCandidates(values: DoubleArray): List<Region> {
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
