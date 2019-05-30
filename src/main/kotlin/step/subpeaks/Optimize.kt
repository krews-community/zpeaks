package step.subpeaks

import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer
import org.apache.commons.math3.linear.RealVector
import org.apache.commons.math3.util.Precision

data class OptimizeResults<T : GaussianParameters>(
    val parameters: List<T>,
    val error: Double,
    val iterations: Int
)

abstract class Optimizer<T : GaussianParameters> {

    fun optimize(values: DoubleArray, candidateGaussians: List<CandidateGaussian<T>>,
                 gaussians: List<T>): OptimizeResults<T> {
        val avg = values.average()

        val optimizer = LevenbergMarquardtOptimizer()
            .withInitialStepBoundFactor(0.05)
            .withCostRelativeTolerance(avg * 1e-4)
            .withParameterRelativeTolerance(avg * 1e-8)
            .withOrthoTolerance(avg * 1e-3)
            .withRankingThreshold(Precision.SAFE_MIN)

        val initialParameters = parametersToArray(gaussians)

        val problem = LeastSquaresBuilder()
            .maxEvaluations(1000)
            .maxIterations(1000)
            .model(
                { params -> calculateCurve(arrayToParameters(params), values.size) },
                { params -> calculateJacobian(params, values.size) }
            )
            .parameterValidator { params -> validateParameters(params, candidateGaussians) }
            .start(initialParameters)
            .target(values)
            .build()

        val optimum = optimizer.optimize(problem)
        val optimizedParameters = arrayToParameters(optimum.point.toArray())
        return OptimizeResults(optimizedParameters, optimum.rms, optimum.iterations)
    }

    fun calculateCurve(parameters: List<T>, curveLength: Int, start: Int = 0): DoubleArray {
        val curve = DoubleArray(curveLength) { 0.0 }
        for (j in 0 until curveLength) {
            for (gaussian in parameters) {
                curve[j] += calculateCurveValue(j+start, gaussian)
            }
        }
        return curve
    }

    protected abstract fun parametersToArray(parameters: List<T>): DoubleArray
    protected abstract fun arrayToParameters(array: DoubleArray): List<T>
    protected abstract fun calculateCurveValue(x: Int, gaussian: T): Double
    protected abstract fun calculateJacobian(parameters: DoubleArray, curveLength: Int): Array<DoubleArray>
    protected abstract fun validateParameters(parameters: RealVector,
                                              candidateGaussians: List<CandidateGaussian<T>>): RealVector

}