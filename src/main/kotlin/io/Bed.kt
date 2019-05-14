package io

import model.Region
import mu.KotlinLogging
import step.Peak
import step.subpeaks.*
import java.nio.ByteBuffer
import java.nio.file.*
import java.util.*

private val log = KotlinLogging.logger {}

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

/**
 * Writes part of bed file for given chromosome and peaks
 */
fun writePeaksBed(path: Path, chr: String, peaks: Iterable<Peak>) {
    log.info { "Writing peaks data for chromosome $chr to bed file $path" }
    Files.newBufferedWriter(path).use { writer ->
        for (peak in peaks) {
            val region = peak.region
            val name = bedPeakName(chr, region)
            writer.write("$chr\t${region.start}\t${region.end}\t$name\t${peak.score}\n")
        }
    }
    log.info { "Chromosome $chr peaks data write complete!" }
}

/**
 * Writes part of bed file for given chromosome and sub-peaks
 */
fun writeSkewSubPeaksBed(path: Path, chr: String, subPeaks: Iterable<SkewSubPeak>) {
    log.info { "Writing skew sub-peaks data for chromosome $chr to bed file $path" }
    Files.newBufferedWriter(path).use { writer ->
        for (subPeak in subPeaks) {
            val region = subPeak.region
            val name = bedPeakName(chr, region)
            val score = "${subPeak.gaussianParameters.amplitude}#${subPeak.gaussianParameters.shape}"
            writer.write("$chr\t${region.start}\t${region.end}\t$name\t$score\n")
        }
    }
    log.info { "Chromosome $chr skew sub-peaks data write complete!" }
}