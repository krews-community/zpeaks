package util

import java.nio.file.*
import kotlin.streams.toList

val TEST_BAM_PATH = getResourcePath("ENCFF375IJW.chr22.bam")

val TEST_BAM_2_PATH = getResourcePath("ERR1299114.bam")

val MULTI_BAM_1_PATH = getResourcePath("ENCFF165CHM.chr22_30m-32m.bam")
val MULTI_BAM_2_PATH = getResourcePath("ENCFF611PME.chr22_30m-32m.bam")
val MULTI_BAM_3_PATH = getResourcePath("ENCFF708IUW.chr22_30m-32m.bam")
val MULTI_BAM_4_PATH = getResourcePath("ENCFF884TCQ.chr22_30m-32m.bam")

val MANY_BAM_DIR_PATH = getResourcePath("ct_bams")
val MANY_BAM_PATHS = Files.walk(MANY_BAM_DIR_PATH)
    .filter { Files.isRegularFile(it) && it.toString().endsWith(".bam") }
    .toList()

const val CHR_22 = "chr22"
const val CHR_22_SIZE = 48_000_000
val CHR_FILTER_PATH = getResourcePath("chr_filter.txt")

val MULTICHROM_PILEUP_BEDGRAPH = getResourcePath("multiChrom.pileup.bedGraph")

interface Junk

fun getResourcePath(relativePath: String): Path {
    val url = Junk::class.java.classLoader.getResource(relativePath)
    return Paths.get(url.toURI())
}

fun Path.filenameWithoutExtension(): String {
    val filename = this.fileName.toString()
    return filename.substring(0, filename.lastIndexOf('.'))
}

fun Path.copyToAndDelete(to: Path): Path {
    Files.copy(this, to, StandardCopyOption.REPLACE_EXISTING)
    Files.deleteIfExists(this)
    return to
}