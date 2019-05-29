package io

import gnu.trove.list.array.*
import model.SignalData
import mu.KotlinLogging
import org.jetbrains.bio.big.*
import java.nio.file.*
import kotlin.math.*


private val log = KotlinLogging.logger {}

data class SignalSection(val chr: String, val start: Int, val span: Int, val value: Number) {
    val end: Int get() = start + span
}

enum class SignalOutputFormat { WIG, BIG_WIG, BED_GRAPH }

fun createSignalFile(signalOut: Path, format: SignalOutputFormat, data: Map<String, SignalData>, signalResolution: Int = 1) {
    when(format) {
        SignalOutputFormat.WIG -> { writeWig(signalOut, data, signalResolution) }
        SignalOutputFormat.BIG_WIG -> { writeBigWig(signalOut, data, signalResolution) }
        SignalOutputFormat.BED_GRAPH -> { writeBedGraph(signalOut, data, signalResolution) }
    }
}

fun writeWig(path: Path, data: Map<String, SignalData>, signalResolution: Int = 1) {
    log.info { "Writing pile-up data to wig file $path" }
    Files.newBufferedWriter(path).use { writer ->
        writer.write("track type=wiggle_0\n")
        var lastSpan: Int? = null
        iteratePileUpSections(data, signalResolution) { section ->
            if (lastSpan != section.span) {
                writer.write("variableStep chrom=${section.chr} span=${section.span}\n")
                lastSpan = section.span
            }
            writer.write("${section.start} ${section.value}\n")
        }
    }
    log.info { "Wig file write complete!" }
}

fun writeBigWig(path: Path, dataByChr: Map<String, SignalData>, signalResolution: Int = 1) {
    log.info { "Writing pile-up data to bigwig file $path" }
    val wigSections = mutableListOf<WigSection>()
    iteratePileUpSections(dataByChr, signalResolution) { section ->
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

fun writeBedGraph(path: Path, dataByChr: Map<String, SignalData>, signalResolution: Int = 1) {
    log.info { "Writing pile-up data to bed-graph file" }
    Files.newBufferedWriter(path).use { writer ->
        writer.write("track type=bedGraph\n")
        iteratePileUpSections(dataByChr, signalResolution) { section ->
            writer.write("${section.chr}\t${section.start}\t${section.end}\t${section.value}\n")
        }
    }
    log.info { "Bed-graph file write complete!" }
}

fun iteratePileUpSections(dataByChr: Map<String, SignalData>, signalResolution: Int, processSection: (section: SignalSection) -> Unit) {
    for ((chr, data) in dataByChr) {
        var currentStart: Int? = null
        var currentValue: Number? = null
        for (index in 0 until data.chrLength) {
            val value = roundValue(data[index], signalResolution)
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

fun roundValue(value: Number, signalResolution: Int): Double {
    return floor(value.toDouble() * (10.0).pow(signalResolution)) / (10.0).pow(signalResolution)
}
