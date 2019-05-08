package step.subpeaks

import model.*
import mu.KotlinLogging
import org.apache.commons.math3.linear.*
import org.apache.commons.math3.linear.MatrixUtils.*
import org.apache.commons.math3.util.FastMath.*
import util.*

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
    val candidates = scaleSpaceZeros(scaledValues)

    /*
     * here we compute guesses for the Gaussians' parameters based on zero crossing locations
     * the means and standard deviations are computed with simple equations
     * the amplitudes are fit to the data using least squares after the means and standard deviations are computed
     */
    val gaussians = candidates.map { initParameters(it) } as MutableList

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

    log.debug { "values: $values" }
    log.debug { "candidate regions: $candidates" }
    log.debug { "m: ${m.map { it.toList() }.toList()}" }
    log.debug { "d: ${d.toList()}" }
    val r = LUDecomposition(createRealMatrix(m)).solver.solve(createRealVector(d))

    for (j in 0 until gaussians.size) {
        val s = if (r.getEntry(j) < 0) -1.0 else 1.0
        // Set Amplitude
        gaussians[j].amplitude = sqrt(abs(r.getEntry(j) * gaussians[j].stdDev) * s)
    }
    // Sort gaussians by amplitude
    gaussians.sortBy { it.amplitude * it.amplitude / it.stdDev }

    /*
     * Perform the actual fit of the curve to a sum of gaussians.size Gaussians
     * optimization is performed with the Levenberg-Marquardt algorithm
     */
    var bestOptimized: OptimizeResults<T>? = null
    for (j in 1 until gaussians.size) {
        val currentGuess = gaussians.subList(0, j)
        val optimizeResults = optimize(scaledValues, currentGuess, 0.1)
        if (optimizeResults.error < 0.05) {
            bestOptimized = optimizeResults
            break
        }
        if (bestOptimized == null || optimizeResults.error < bestOptimized.error) {
            bestOptimized = optimizeResults
        }
    }

    val sqrtAvg = sqrt(avg)
    // Bring parameters back to original scale
    bestOptimized!!.parameters.forEach { it.amplitude *= sqrtAvg }

    return bestOptimized.parameters.map {
        SubPeak(parametersToRegion(it, start), parametersToScore(it), it)
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
fun scaleSpaceZeros(values: DoubleArray): List<Region> {
    val candidates = mutableListOf<Region>()

    // Start with a kernel near the square root of the vector size and loop down to 2
    val kernel = Math.floor(sqrt(values.size.toDouble()) * 1.1).toInt()
    for (kernelSize in kernel downTo 2 step 1) {

        // Smooth at this kernel width
        val blurred = scaleSpaceSmooth(values, kernelSize.toDouble())
        log.info { "blurred: $blurred" }

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
