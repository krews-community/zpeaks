package io

import org.jetbrains.bio.big.WigPrinter
import step.PileUp
import java.nio.file.Files
import java.nio.file.Path

// TODO
/*
fun writeWig(path: Path, name: String, pileUps: Map<String, PileUp>) {
    WigPrinter(Files.newBufferedWriter(path), name).use {
        pileUps.forEach { (chr, pileUp) ->
            for (chrIndex in 0 until pileUp.chromosomeLength) {
                val pileUpValue = pileUp.values[chrIndex] ?: 0
            }
        }
    }
}*

 */