package io

import model.Region
import mu.KotlinLogging
import step.Peak
import step.subpeaks.*
import java.nio.ByteBuffer
import java.nio.file.*
import java.util.*

private val log = KotlinLogging.logger {}

private val Double.printValue: String get() = "%.2f".format(this)

/**
 * Writes part of bed file for given chromosome and peaks
 */
fun writePeaksBed(path: Path, peaks: Map<String, Iterable<Peak>>) {
    log.info { "Writing peaks data to bed file $path" }
    Files.newBufferedWriter(path).use { writer ->
        for ((chr, chrPeaks) in peaks) {
            for (peak in chrPeaks) {
                val region = peak.region
                val name = bedPeakName(chr, region)
                writer.write("$chr\t${region.start}\t${region.end}\t$name\t${peak.score.printValue}\n")
            }
        }
    }
    log.info { "Peaks data write complete!" }
}

/**
 * Writes part of bed file for given chromosome and sub-peaks
 */
fun writeSkewSubPeaksBed(path: Path, subPeaks: Map<String, Iterable<SubPeak<SkewGaussianParameters>>>) {
    log.info { "Writing skew sub-peaks data to bed file $path" }
    Files.newBufferedWriter(path).use { writer ->
        for ((chr, chrSubPeaks) in subPeaks) {
            for (subPeak in chrSubPeaks) {
                val region = subPeak.region
                val name = bedPeakName(chr, region)
                val amplitude = subPeak.gaussianParameters.amplitude
                val shape = subPeak.gaussianParameters.shape
                val score = "${amplitude.printValue}#${shape.printValue}"
                writer.write("$chr\t${region.start}\t${region.end}\t$name\t$score\n")
            }
        }
    }
    log.info { "Skew sub-peaks data write complete!" }
}

fun writeStandardSubPeaksBed(path: Path, subPeaks: Map<String, Iterable<SubPeak<StandardGaussianParameters>>>) {
    log.info { "Writing Standard Sub-Peaks data to bed file $path" }
    Files.newBufferedWriter(path).use { writer ->
        for ((chr, chrSubPeaks) in subPeaks) {
            for (subPeak in chrSubPeaks) {
                val region = subPeak.region
                val name = bedPeakName(chr, region)
                val amplitude = subPeak.gaussianParameters.amplitude
                val score = amplitude.printValue
                writer.write("$chr\t${region.start}\t${region.end}\t$name\t$score\n")
            }
        }
    }
    log.info { "Standard Sub-Peaks data write complete!" }
}

/**
 * Creates a small-ish name unique for each chromosome range
 * The name consists the chromosome name plus a base64 encoded (and therefor shrunken) version of the range
 */
fun bedPeakName(chr: String, region: Region): String {
    val base64Range = Base64.getEncoder()
        .encodeToString(region.start.toByteArray() + region.end.toByteArray())
        .trim('=')
    return chr + base64Range
}

private fun Int.toByteArray(): ByteArray = ByteBuffer.allocate(4).putInt(this).array()