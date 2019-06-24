package model

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