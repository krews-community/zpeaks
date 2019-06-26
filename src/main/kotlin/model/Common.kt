package model

import io.SignalOutputFormat
import step.PileUpOptions
import java.nio.file.Path

enum class Strand {
    PLUS, MINUS, BOTH
}

data class Region(val start: Int, val end: Int)

interface SignalData {
    val chrLength: Int
    operator fun get(bp: Int): Number
}

data class PileUpInput(val bam: Path, val options: PileUpOptions)

enum class FitMode { SKEW, STANDARD }
data class SignalOutput(
    val path: Path,
    val type: SignalOutputType,
    val format: SignalOutputFormat,
    val signalResolution: Int = 1
)
enum class SignalOutputType { RAW, SMOOTHED }