import htsjdk.samtools.SamReaderFactory
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test
import step.prepBams
import util.*

class BamPrepTests {

    @Test
    fun `Test Prep Bams`() {
        val humanChrs = ((1 until 22).map { it.toString() } + "X" + "Y" + "M").map { "chr$it" }
        val chrLengths = prepBams(MANY_BAM_PATHS, humanChrs)
        assertThat(chrLengths.keys.toSet()).isEqualTo(humanChrs.toSet())

        var sum = 0
        SamReaderFactory.make().open(MANY_BAM_PATHS[0]).use { reader ->
            val i = reader.indexing().index
            for(record in reader.query("chr22", 0, 0, false)) {
                sum++
            }
        }
        assertThat(sum).isGreaterThan(0)
    }

}