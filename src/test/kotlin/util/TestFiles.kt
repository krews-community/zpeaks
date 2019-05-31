package util

import java.nio.file.*

val TEST_BAM_PATH = getResourcePath("ENCFF375IJW.chr22.bam")
const val TEST_BAM_CHR = "chr22"

val TEST_BAM_2_PATH = getResourcePath("ERR1299114.bam")

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