package io

import java.nio.file.*

fun parseChrFilter(path: Path): Map<String, IntRange?> {
    val chrFilter = mutableMapOf<String, IntRange?>()
    Files.newBufferedReader(path).forEachLine { line ->
        try {
            val splitLine = line.split("\\s+".toRegex())
            val chr = splitLine[0]
            if (splitLine.size == 1) {
                chrFilter[chr] = null
                return@forEachLine
            }
            val rangeParts = splitLine[1].split("-")
            val range = rangeParts[0].toInt() until rangeParts[1].toInt()
            chrFilter[chr] = range
        } catch (e: Exception) {
            throw Exception("Invalid chromosome filter file entry on line:\n$line", e)
        }
    }
    return chrFilter
}