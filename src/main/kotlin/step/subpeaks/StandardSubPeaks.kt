package step.subpeaks

import model.*
import org.apache.commons.math3.linear.RealVector
import kotlin.math.*


data class StandardGaussianParameters(
    override var amplitude: Double,
    override var mean: Double,
    override var stdDev: Double
) : GaussianParameters

object StandardFitter : Fitter<StandardGaussianParameters>("Fit Standard Sub-Peaks", StandardOptimizer) {
    override fun initParameters(region: Region) = StandardGaussianParameters(
        amplitude = 0.0,
        mean = (region.start + region.end) / 2.0,
        stdDev = (region.end - region.start) / 2.0
    )

    override fun parametersToRegion(parameters: StandardGaussianParameters) = Region(
        start = (parameters.mean - parameters.stdDev).toInt(),
        end = (parameters.mean + parameters.stdDev).toInt()
    )
}

object StandardOptimizer : Optimizer<StandardGaussianParameters>() {

    override fun parametersToArray(parameters: List<StandardGaussianParameters>): DoubleArray {
        val paramArray = DoubleArray(parameters.size * 3)
        for ((i, gaussian) in parameters.withIndex()) {
            paramArray[i * 3] = gaussian.amplitude
            paramArray[i * 3 + 1] = gaussian.mean
            paramArray[i * 3 + 2] = gaussian.stdDev
        }
        return paramArray
    }

    override fun arrayToParameters(array: DoubleArray): List<StandardGaussianParameters> {
        val parameters = mutableListOf<StandardGaussianParameters>()
        for (j in 0 until array.size step 3) {
            parameters += StandardGaussianParameters(array[j], array[j + 1], array[j + 2])
        }
        return parameters
    }

    override fun calculateCurveValue(x: Int, gaussian: StandardGaussianParameters) = with(gaussian) {
        amplitude * amplitude / stdDev * exp(-(x - mean) * (x - mean) / stdDev / stdDev / 2)
    }

    override fun calculateJacobian(parameters: DoubleArray, curveLength: Int): Array<DoubleArray> {
        val jacobian = Array(curveLength) { DoubleArray(parameters.size) }
        for (j in 0 until curveLength) {
            for (k in 0 until parameters.size step 3) {
                val a = parameters[k] // amplitude
                val m = parameters[k + 1] // mean
                val u = parameters[k + 2] // standard deviation
                val exp = exp(-(j - m) * (j - m) / (2.0 * u * u))
                // partial with respect to a for this Gaussian
                jacobian[j][k] = 2 * a * exp / u
                // partial with respect to m
                jacobian[j][k+1] = a * a * (j - m) * exp / (u * u * u)
                // partial with respect to u
                jacobian[j][k+2] = a * a * exp * (j * j - 2 * m * j - u * u + m * m) / (u * u * u * u)
            }
        }
        return jacobian
    }

    private const val MAX_STD_DEV = 400

    override fun validateParameters(
        parameters: RealVector,
        candidateGaussians: List<CandidateGaussian<StandardGaussianParameters>>
    ): RealVector {
        val validated = parameters.copy()
        for (j in 0 until candidateGaussians.size) {
            val amplitude = parameters.getEntry(j * 3)
            val mean = parameters.getEntry(j * 3 + 1)
            val stdDev = parameters.getEntry(j * 3 + 2)

            if (amplitude <= 0) {
                validated.setEntry(j * 3, candidateGaussians[j].parameters.amplitude)
            }

            if (mean < candidateGaussians[j].region.start || mean > candidateGaussians[j].region.end) {
                validated.setEntry(j * 3 + 1, candidateGaussians[j].parameters.mean)
            }
            if (stdDev <= 0 || stdDev > MAX_STD_DEV) {
                validated.setEntry(j * 3 + 2, candidateGaussians[j].parameters.stdDev)
            }
        }
        return validated
    }

}