package step

import model.ChromBounds
import model.SignalData
import java.nio.file.Path
import kotlin.math.sqrt

class PileUp(
        // Pile-up values in chromosome.
        private val values: FloatArray,
        // Chromosome as named in alignment file
        override val chr: String,
        // Chromosome length pulled directly from BAM File
        val chrLength: Int,
        // Chromosome range to be used downstream
        override val range: IntRange,
        // Sum of all pile-up values in chromosome calculated on the fly and cached for efficiency
        val sum: Double
) : SignalData {
    override operator fun get(bp: Int): Float = values[bp]

    private val scalingFactor: Float = if(sum == 0.0) 1.0F else sqrt(sum).toFloat()
    fun scaledValue(bp: Int): Float = this[bp] / scalingFactor
}

interface PileUpRunner {
    fun getChromsWithBounds(chrFilter: Map<String, IntRange?>? = null): Map<String, ChromBounds>

    fun getPileUps(chr: String, chrLength: Int, range: IntRange): Sequence<PileUp>
}
