import org.assertj.core.api.Assertions
import org.junit.jupiter.api.Test
import step.BedGraphPileUpRunner
import util.MULTICHROM_PILEUP_BEDGRAPH
import util.isBedGraph

class BedGraphPileUpRunnerTests {

    @Test
    fun `Test Chrom Offsets`() {
        val runner = BedGraphPileUpRunner(listOf(MULTICHROM_PILEUP_BEDGRAPH))
        val chromOffsets = runner.chromOffsets
        val expected = listOf<Triple<String, Int, Long>>(
            Triple("chr20", 64332156, 20),
            Triple("chr21", 46700008, 21417232),
            Triple("chr22", 50808468, 29776620)
        )
        // Test in reverse since later offsets will cause earlier sizes to be wrong
        for (chr in expected.reversed()) {
            Assertions.assertThat(chromOffsets).containsKey(chr.first)
            val offsets = chromOffsets[chr.first]!!
            Assertions.assertThat(offsets).containsKey(MULTICHROM_PILEUP_BEDGRAPH)
            val offsetsBedGraph = offsets[MULTICHROM_PILEUP_BEDGRAPH]
            Assertions.assertThat(offsetsBedGraph).isEqualTo(Pair(chr.second, chr.third))
        }
    }

    @Test
    fun `Test getPileUps`() {
        val runner = BedGraphPileUpRunner(listOf(MULTICHROM_PILEUP_BEDGRAPH))
        val chrsWithBounds = runner.getChromsWithBounds()
        val pileUps = runner.getPileUps("chr22", chrsWithBounds["chr22"]!!.length, SAMPLE_RANGE_LARGE).toList()
        Assertions.assertThat(pileUps.size).isEqualTo(1)
    }

    @Test
    fun `Test isBedGraph`() {
        Assertions.assertThat(isBedGraph(MULTICHROM_PILEUP_BEDGRAPH)).isTrue()
    }
}