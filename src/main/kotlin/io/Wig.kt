package io

import gnu.trove.list.array.*
import model.SignalData
import mu.KotlinLogging
import org.jetbrains.bio.big.*
import java.nio.file.*


private val log = KotlinLogging.logger {}

data class SignalSection(val chr: String, val start: Int, val span: Int, val value: Number) {
    val end: Int get() = start + span
}

enum class SignalOutputFormat { WIG, BIG_WIG, BED_GRAPH }

fun createSignalFile(signalOut: Path, format: SignalOutputFormat, data: Map<String, SignalData>) {
    when(format) {
        SignalOutputFormat.WIG -> { writeWig(signalOut, data) }
        SignalOutputFormat.BIG_WIG -> { writeBigWig(signalOut, data) }
        SignalOutputFormat.BED_GRAPH -> { writeBedGraph(signalOut, data) }
    }
}

fun writeWig(path: Path, data: Map<String, SignalData>) {
    log.info { "Writing pile-up data to wig file $path" }
    Files.newBufferedWriter(path).use { writer ->
        writer.write("track type=wiggle_0\n")
        var lastSpan: Int? = null
        iteratePileUpSections(data) { section ->
            if (lastSpan != section.span) {
                writer.write("variableStep chrom=${section.chr} span=${section.span}\n")
                lastSpan = section.span
            }
            writer.write("${section.start} ${section.value}\n")
        }
    }
    log.info { "Wig file write complete!" }
}

fun writeBigWig(path: Path, dataByChr: Map<String, SignalData>) {
    log.info { "Writing pile-up data to bigwig file $path" }
    val wigSections = mutableListOf<WigSection>()
    iteratePileUpSections(dataByChr) { section ->
        wigSections += VariableStepSection(
            chrom = section.chr,
            span = section.span,
            positions = TIntArrayList(intArrayOf(section.start)),
            values = TFloatArrayList(floatArrayOf(section.value.toFloat()))
        )
    }
    val chromSizes = dataByChr.map { (chr, data) -> chr to data.chrLength }
    BigWigFile.write(wigSections, chromSizes, path)
    log.info { "Bigwig file write complete!" }
}

fun writeBedGraph(path: Path, dataByChr: Map<String, SignalData>) {
    log.info { "Writing pile-up data to bed-graph file" }
    Files.newBufferedWriter(path).use { writer ->
        writer.write("track type=bedGraph\n")
        iteratePileUpSections(dataByChr) { section ->
            writer.write("${section.chr}\t${section.start}\t${section.end}\t${section.value}\n")
        }
    }
    log.info { "Bed-graph file write complete!" }
}

fun iteratePileUpSections(dataByChr: Map<String, SignalData>, processSection: (section: SignalSection) -> Unit) {
    for ((chr, data) in dataByChr) {
        var currentStart: Int? = null
        var currentValue: Number? = null
        for (index in 0 until data.chrLength) {
            val value = data[index]
            if (currentValue != value) {
                if (currentValue != null) {
                    processSection(SignalSection(chr, currentStart!!, index - currentStart, currentValue))
                }
                if (value.toDouble() > 0.0) {
                    currentStart = index
                    currentValue = value
                } else {
                    currentStart = null
                    currentValue = null
                }
            }
        }
    }
}