package model

enum class Strand {
    PLUS, MINUS, BOTH
}

data class Region(val start: Int, val end: Int)

interface SignalData {
    val chrLength: Int
    operator fun get(bp: Int): Number
}