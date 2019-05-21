package step.subpeaks

import model.Region
import org.apache.commons.math3.linear.LUDecomposition
import org.apache.commons.math3.linear.MatrixUtils.*
import org.apache.commons.math3.util.FastMath.*
import util.gaussianDistribution
import util.scaleSpaceSmooth


fun <T : GaussianParameters> findCandidates(
    data: List<Double>,
    initParameters: (region: Region) -> T
): List<CandidateGaussian<T>> {
    var zeroCrossings = mutableListOf<Int>()

    // Smooth at this kernel width
    val blurred = scaleSpaceSmooth(data, 1.0)

    // Get the zero crossings for this kernel size
    var cSign = blurred[0] > 0
    for (j in 1 until blurred.size) {
        if (cSign != blurred[j] > 0) {
            cSign = blurred[j] > 0
            zeroCrossings.add(j)
        }
    }

    if (zeroCrossings.size % 2 == 1) {
        val first = data.first()
        val last = data.last()
        if (first > last) {
            zeroCrossings = zeroCrossings.subList(1, zeroCrossings.size)
        }
        if (last > first) {
            zeroCrossings = zeroCrossings.subList(0, zeroCrossings.size-1)
        }
    }

    var candidateRegions = zeroCrossingsToRegions(zeroCrossings)
    var candidateGaussians = candidateGaussians(data, candidateRegions, initParameters)
    if (candidateGaussians.any { it.parameters.amplitude < 0 }) {
        zeroCrossings = zeroCrossings.subList(1, zeroCrossings.size-1)
        candidateRegions = zeroCrossingsToRegions(zeroCrossings)
        candidateGaussians = candidateGaussians(data, candidateRegions, initParameters)
    }
    return candidateGaussians
}

private fun zeroCrossingsToRegions(zeroCrossings: List<Int>): List<Region> {
    val candidates = mutableListOf<Region>()
    for (i in 0 until zeroCrossings.size step 2) {
        candidates.add(Region(zeroCrossings[i], zeroCrossings[i+1]))
    }
    return candidates
}

fun <T : GaussianParameters> candidateGaussians(
    values: List<Double>,
    candidateRegions: List<Region>,
    initParameters: (region: Region) -> T
): List<CandidateGaussian<T>> {
    /*
     * here we compute guesses for the Gaussians' parameters based on zero crossing locations
     * the means and standard deviations are computed with simple equations
     * the amplitudes are fit to the data using least squares after the means and standard deviations are computed
     */
    val candidateGaussians = candidateRegions
        .map { CandidateGaussian(it, initParameters(it)) } as MutableList

    /*
     * Solve the system of linear equations to get the amplitude guesses
     */
    val lim = values.size
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
            value += distributions[j][n] * values[n]
        }
        d[j] = value
    }

    val r = LUDecomposition(createRealMatrix(m)).solver.solve(createRealVector(d))

    for (j in 0 until candidateGaussians.size) {
        val s = if (r.getEntry(j) < 0) -1.0 else 1.0
        // Set Amplitude
        candidateGaussians[j].parameters.amplitude =
            sqrt(abs(r.getEntry(j)) * candidateGaussians[j].parameters.stdDev) * s
    }

    return candidateGaussians
}