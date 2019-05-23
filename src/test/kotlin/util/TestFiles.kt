package util

import java.nio.file.*

val TEST_BAM_PATH = getResourcePath("ENCFF375IJW.chr22.bam")
const val TEST_BAM_CHR = "chr22"

interface Junk

fun getResourcePath(relativePath: String): Path {
    val url = Junk::class.java.classLoader.getResource(relativePath)
    return Paths.get(url.toURI())
}