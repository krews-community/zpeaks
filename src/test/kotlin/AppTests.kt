import htsjdk.samtools.SamReaderFactory
import org.junit.jupiter.api.Test
import util.getResourcePath

class AppTests {

    @Test
    fun x() {
        val testBamPath = getResourcePath("test.bam")
        val bamReader = SamReaderFactory.make().open(testBamPath)
        val sdict = bamReader.fileHeader.sequenceDictionary
        val chromosomeLengths: Map<String, Int> = bamReader.fileHeader.sequenceDictionary.sequences
            .map { it.sequenceName to it.sequenceLength }.toMap()
        bamReader.use { reader ->
            reader.forEach { record ->
                val x = record.referenceName
            }
        }
    }
}