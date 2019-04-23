package util

import java.lang.Math.*
import kotlin.math.pow

fun gaussianDistribution(a: Double, x: Double, u: Double, length: Int): List<Double> {
    return (0..length).map { i -> a * exp(-(i - x) * (i - x) / (u * u) / 2.0) }
}

/**
 * Smooths the input vector using the second derivative of a Gaussian distribution
 * of the given width.
 */
fun scaleSpaceSmooth(input: List<Double>, kernelWidth: Double): List<Double> {
    val gaussianKernel = sampleSDev(kernelWidth, input.size)
    return convolve(input, gaussianKernel)
}

/**
 * Perform a 1d convolution with output list same length as input
 */
fun convolve(input: List<Double>, kernel: List<Double>): List<Double> {
    val output = mutableListOf<Double>()
    val kernelCenter = kernel.size / 2
    for (inputIndex in 0 until input.size) {
        var outputVal = 0.0
        for ((kernelIndex, kernelVal) in kernel.withIndex()) {
            val neighbor = inputIndex + kernelIndex - kernelCenter
            if (neighbor >= 0 && neighbor < input.size) {
                outputVal += input[neighbor] * kernelVal
            }
        }
        output[inputIndex] = outputVal
    }
    return output
}

/**
 * Samples the second derivative of a Gaussian distribution at
 * a given number of points.
 *
 * stdev: the standard deviation of the distribution.
 * size: the size of the output vector generated will be size * 2 + 1,
 *       with the mean at the center.
 */
fun sampleSDev(stdDev: Double, size: Int): List<Double> {
    val a = 1.0 / stdDev / SQRT2PI
    val output = mutableListOf<Double>()
    for (i in 0 until size) {
        val value = a *
                exp(-i.pow(2) / stdDev.pow(2)) *
                (i.pow(2) - stdDev.pow(2)) /
                stdDev.pow(4)
        output[size + i] = value
        output[size - i] = value
    }
    return output
}