package io

import model.Region
import step.Peak
import step.subpeaks.*
import java.nio.file.*
import java.util.*

fun bedPeakName(chr: String, region: Region): String {
    return chr + Base64.getEncoder().encode(byteArrayOf(region.start.toByte(), region.end.toByte()))
}

/**
 * Writes part of bed file for given chromosome and peaks
 */
fun writePeaksBed(path: Path, chr: String, peaks: Iterable<Peak>) {
    Files.newBufferedWriter(path).use { writer ->
        for (peak in peaks) {
            val region = peak.region
            val name = bedPeakName(chr, region)
            writer.write("$chr\t${region.start}\t${region.end}\t$name\t${peak.score}\n")
        }
    }
}

/**
 * Writes part of bed file for given chromosome and sub-peaks
 */
fun writeSkewSubPeaksBed(path: Path, chr: String, subPeaks: Iterable<SkewSubPeak>) {
    Files.newBufferedWriter(path).use { writer ->
        for (subPeak in subPeaks) {
            val region = subPeak.region
            val name = bedPeakName(chr, region)
            val score = "${subPeak.gaussianParameters.amplitude}#${subPeak.gaussianParameters.shape}"
            writer.write("$chr\t${region.start}\t${region.end}\t$name\t$score\n")
        }
    }
}