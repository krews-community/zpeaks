package step

import com.google.common.util.concurrent.AtomicDoubleArray
import model.*
import mu.KotlinLogging
import util.*
import java.util.stream.IntStream
import kotlin.random.Random

private val log = KotlinLogging.logger {}

const val BACKGROUND_LIMIT = 1000

data class Background(val average: Double, val stdDev: Double)

class PDF(
    private val values: AtomicDoubleArray,
    val background: Background,
    override val chrLength: Int
): SignalData {
    override operator fun get(bp: Int): Double = values[bp]
}

fun runSmooth(pileUps: MutableMap<String, PileUp>, smoothing: Double, normalizePDF: Boolean): Map<String, PDF> {
    log.info { "Performing smoothing of raw pile up data for ${pileUps.size} chromosomes..." }
    val pdfs = mutableMapOf<String, PDF>()
    for ((chr, pileUp) in pileUps) {

        log.info { "Calculating PDF for chromosome $chr pileup data..." }
        val pdf = pdf(chr, pileUp, smoothing, normalizePDF)
        log.info { "Chromosome $chr PDF completed with background ${pdf.background}" }

        if (0.0 == pdf.background.average || 0.0 == pdf.background.stdDev) {
            log.warn { "One or more background parameters for chromosome $chr was zero. skipping..." }
            continue
        }
        pdfs[chr] = pdf

        // Remove pile-up data so it can be garbage collected and memory can be reclaimed.
        pileUps.remove(chr)
    }
    log.info { "Smoothing complete!" }
    return pdfs
}

/**
 * Using FSeq algorithm, calculate the probability density function for the given pile-up data
 */
fun pdf(chr: String, pileUp: PileUp, bandwidth: Double, normalizePDF: Boolean, onRange: IntRange? = null): PDF {
    val windowSize = windowSize(bandwidth)
    val lookupTable = lookupTable(normalizePDF, windowSize, bandwidth, pileUp.sum)
    val pdfValues = AtomicDoubleArray(pileUp.chrLength)
    logProgress("Creating PDF for $chr", pileUp.chrLength) { tracker ->
        val start = onRange?.start ?: 0
        val end = onRange?.endInclusive ?: pileUp.chrLength
        IntStream.range(start, end).parallel().forEach { chrIndex ->
            val pileUpValue = pileUp[chrIndex]
            if (pileUpValue == 0) return@forEach
            pdfValues.addAndGet(chrIndex, pileUpValue * lookupTable[0])
            for (i in 1 until windowSize) {
                if (chrIndex + i < pileUp.chrLength) {
                    pdfValues.addAndGet(chrIndex + i, pileUpValue * lookupTable[i])
                }
                if (chrIndex - i > 0) {
                    pdfValues.addAndGet(chrIndex - i, pileUpValue * lookupTable[i])
                }
            }
            tracker.incrementAndGet()
        }
    }

    val background = background(lookupTable, pileUp.sum, windowSize, pileUp.chrLength)
    return PDF(pdfValues, background, pileUp.chrLength)
}

private fun windowSize(bandwidth: Double): Int {
    return (Math.sqrt(Math.log(Double.MIN_VALUE * SQRT2PI) * -2.0) * bandwidth).toInt()
}

/**
 * Compute the background. We get the average number of regions per window, then get n
 * random distributions of this many reads across the window. We repeat 1000 times; for
 * each repeat, we compute the PDF at the center; at high numbers of repeats, the
 * distribution of PDFs should be near normal.
 */
private fun background(lookupTable: List<Double>, numBins: Int, windowSize: Int, chrLength: Int): Background {
    val averageN = numBins * windowSize / chrLength
    val backgroundDist = mutableListOf<Double>()
    for (i in 0 until BACKGROUND_LIMIT) {
        var bgVal = 0.0
        for (j in 0 until averageN) {
            val x = Random.nextInt(0, windowSize / 2)
            bgVal += lookupTable[x]
        }
        backgroundDist += bgVal
    }

    // Compute step.background average and standard deviation from the distribution.
    val average = backgroundDist.average()
    val stdDev = Math.sqrt(backgroundDist.sumByDouble { Math.pow(it - average, 2.0) } / BACKGROUND_LIMIT)
    return Background(average, stdDev)
}

private fun lookupTable(normalizePDF: Boolean, windowSize: Int, bandwidth: Double, n: Int): List<Double> {
    val lookupTable = if (windowSize > 0) {
        val a = if (!normalizePDF) 1.0 else 1.0 / bandwidth / SQRT2PI
        gaussianDistribution(a, 0.0, bandwidth, windowSize)
    } else {
        listOf(1.0)
    }
    return if (normalizePDF) lookupTable.map { it / n } else lookupTable
}