package io

import gnu.trove.list.array.*
import mu.KotlinLogging
import org.jetbrains.bio.big.*
import step.PileUp
import java.nio.file.*


private val log = KotlinLogging.logger {}

data class PileUpSection(val chr: String, val start: Int, val span: Int, val value: Int) {
    val end: Int get() = start + span
}

fun iteratePileUpSections(pileUps: Map<String, PileUp>, processSection: (section: PileUpSection) -> Unit) {
    for ((chr, pileUp) in pileUps) {
        var currentStart: Int? = null
        var currentValue: Int? = null
        for (chrIndex in 0 until pileUp.chromosomeLength) {
            val pileUpValue = pileUp[chrIndex]
            if (currentValue != pileUpValue) {
                if (currentValue != null) {
                    processSection(PileUpSection(chr, currentStart!!, chrIndex - currentStart, currentValue))
                }
                if (pileUpValue > 0) {
                    currentStart = chrIndex
                    currentValue = pileUpValue
                } else {
                    currentStart = null
                    currentValue = null
                }
            }
        }
    }
}

fun writeWig(path: Path, pileUps: Map<String, PileUp>) {
    log.info { "Writing pile-up data to wig file $path" }
    Files.newBufferedWriter(path).use { writer ->
        writer.write("track type=wiggle_0\n")
        var lastSpan: Int? = null
        iteratePileUpSections(pileUps) { section ->
            if (lastSpan != section.span) {
                writer.write("variableStep chrom=${section.chr} span=${section.span}\n")
                lastSpan = section.span
            }
            writer.write("${section.start} ${section.value}\n")
        }
    }
    log.info { "Wig file write complete!" }
}

fun writeBigWig(path: Path, pileUps: Map<String, PileUp>) {
    log.info { "Writing pile-up data to bigwig file $path" }
    val wigSections = mutableListOf<WigSection>()
    iteratePileUpSections(pileUps) { section ->
        wigSections += VariableStepSection(
            chrom = section.chr,
            span = section.span,
            positions = TIntArrayList(intArrayOf(section.start)),
            values = TFloatArrayList(floatArrayOf(section.value.toFloat()))
        )
    }
    val chromSizes = pileUps.map { (chr, pileUp) -> chr to pileUp.chromosomeLength }
    BigWigFile.write(wigSections, chromSizes, path)
    log.info { "Bigwig file write complete!" }
}

fun writeBedGraph(path: Path, pileUps: Map<String, PileUp>) {
    log.info { "Writing pile-up data to bed-graph file" }
    Files.newBufferedWriter(path).use { writer ->
        writer.write("track type=bedGraph\n")
        iteratePileUpSections(pileUps) { section ->
            writer.write("${section.chr}\t${section.start}\t${section.end}\t${section.value}\n")
        }
    }
    log.info { "Bed-graph file write complete!" }
}