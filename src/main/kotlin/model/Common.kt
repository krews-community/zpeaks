package model

enum class Strand {
    PLUS, MINUS, BOTH
}

data class Region(val start: Int, val end: Int)