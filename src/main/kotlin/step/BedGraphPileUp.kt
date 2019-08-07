package step

import model.*
import mu.KotlinLogging
import util.logProgress
import java.io.RandomAccessFile
import java.nio.file.Path


private val log = KotlinLogging.logger {}

data class BedGraphPileUpRunner(
    val bedGraphs: List<Path>
) : PileUpRunner {
    // pair is size, offset
    val chromOffsets: Map<String, Map<Path, Pair<Int, Long>>> by lazy {
        log.info { "Indexing chromosome offsets in pileup bedgraphs." }
        // We're doing this for two reasons:
        // 1) We need to all the chromosomes and their sizes
        // 2) We want to index the boundaries in the file for faster access
        //
        // This approach uses a binary tree-like approach to find the file offsets between
        val ret = mutableMapOf<String, MutableMap<Path, Pair<Int, Long>>>()
        for (bedGraph in bedGraphs) {
            val file = RandomAccessFile(bedGraph.toFile(), "r")
            var headerEnd = file.filePointer
            while (true) {
                val r = file.readLine()
                if (!r.startsWith("track")) {
                    break
                }
                headerEnd = file.filePointer
            }
            val fileEnd = file.length()
            val offsets = mutableMapOf<String, Long>()
            getOffsets(file, headerEnd, fileEnd, offsets)
            val sizes = getSizes(file, offsets)
            for (chr in offsets.keys) {
                val offset = offsets[chr]!!
                val size = sizes[chr]!!
                ret.getOrPut(chr) { mutableMapOf() }.putIfAbsent(bedGraph, Pair(size, offset))
            }
        }
        ret
    }

    fun getOffsets(file: RandomAccessFile, offsetStart: Long, offsetEnd: Long, offsets: MutableMap<String, Long>, passedLeft: String? = null, passedRight: String? = null) {
        // This wouldn't be a valid file anyways (maybe empty), but this prevents errors/unexpected behavior
        if (offsetEnd - offsetStart < 2) {
            return
        }
        val left = if (passedLeft != null) {
            passedLeft
        } else {
            file.seek(offsetStart)
            val left = file.readLine().split("\t")[0]
            offsets[left] = offsetStart
            left
        }
        val right = if (passedRight != null) {
            passedRight
        } else {
            file.seek(offsetEnd)
            while (true) {
                val nextOffset = file.filePointer - 2
                file.seek(nextOffset)
                if (nextOffset <= offsetStart) {
                    return
                }
                if (file.readByte() == '\n'.toByte()) {
                    break
                }
            }
            val line = file.readLine()
            line.split("\t")[0]
        }
        if (left == right) {
            return
        }
        file.seek((offsetEnd - offsetStart) / 2 + offsetStart)
        while (true) {
            val nextOffset = file.filePointer - 2
            file.seek(nextOffset)
            if (nextOffset <= offsetStart) {
                break
            }
            if (file.readByte() == '\n'.toByte()) {
                break
            }
        }
        val mid = file.filePointer
        if (mid == offsetStart) {
            // The current section is too small to subdivide
            // Just read the entire section
            file.seek(offsetStart)
            var currentOffset = offsetStart
            var lastchrom: String? = null
            while (true) {
                currentOffset = file.filePointer
                if (currentOffset >= offsetEnd) {
                    return
                }
                val r = file.readLine()?.split("\t") ?: return
                if (lastchrom != r[0] && lastchrom != null) {
                    offsets[r[0]] = currentOffset
                }
                lastchrom = r[0]
            }
        }
        while (true) {
            val nextOffset = file.filePointer - 2
            file.seek(nextOffset)
            if (nextOffset <= offsetStart) {
                break
            }
            if (file.readByte() == '\n'.toByte()) {
                break
            }
        }
        val midright = file.readLine().split("\t")[0]
        assert(mid == file.filePointer)
        val midleft = file.readLine().split("\t")[0]
        if (midright != midleft) {
            offsets[midleft] = mid
        }
        getOffsets(file, offsetStart, mid, offsets, left, midright)
        getOffsets(file, mid, offsetEnd, offsets, midleft, right)
    }

    fun getSizes(file: RandomAccessFile, offsets: Map<String, Long>): Map<String, Int> {
        val sizes = mutableMapOf<String, Int>()
        var lastchrom: String? = null
        offsets.entries.sortedBy { it.value }.forEach {
            if (lastchrom == null) {
                lastchrom = it.key
                return@forEach
            }
            file.seek(it.value)
            while (true) {
                val nextOffset = file.filePointer - 2
                file.seek(nextOffset)
                if (file.readByte() == '\n'.toByte()) {
                    break
                }
            }
            sizes[lastchrom!!] = file.readLine().split("\t")[2].toInt()
            lastchrom = it.key
        }
        file.seek(file.length())
        while (true) {
            val nextOffset = file.filePointer - 2
            file.seek(nextOffset)

            if (file.readByte() == '\n'.toByte()) {
                break
            }
        }
        sizes[lastchrom!!] = file.readLine().split("\t")[2].toInt()
        return sizes
    }

    override fun getChromsWithBounds(chrFilter: Map<String, IntRange?>?): Map<String, ChromBounds> {
        return chromOffsets.entries
            .filter { chrFilter == null || chrFilter.contains(it.key) }
            .map {
                val chr = it.key
                val map = it.value
                var length = 0
                for ((_, pair) in map) {
                    val (size, _) = pair
                    length = maxOf(length, size)
                }
                val range =
                    if (chrFilter?.get(chr) != null) chrFilter.getValue(chr)!!
                    else (0 until length)
                chr to ChromBounds(chr, length, range)
            }.toMap()
    }

    override fun getPileUps(chr: String, chrLength: Int, range: IntRange): Sequence<PileUp> = sequence {
        for (bedGraph in bedGraphs) {
            log.info { "Reading pile-up for chromosome $chr on file $bedGraph..." }
            val values = FloatArray(chrLength) { 0.0F }
            var sum = 0.0
            val file = RandomAccessFile(bedGraph.toFile(), "r")
            file.seek(chromOffsets.getValue(chr).getValue(bedGraph).second)
            logProgress("BedGraphPileupReader", chrLength) {tracker ->
                while (true) {
                    val line = file.readLine() ?: break
                    val split = line.split("\t")
                    val chrom = split[0]
                    if (chrom != chr) break
                    val start = split[1].toInt()
                    val end = split[2].toInt()
                    val value = split[3].toFloat()
                    val length = end - start
                    for (base in start until end) {
                        values[base] += value
                    }
                    sum += value * length

                    tracker.set(end)
                }
            }
            yield(PileUp(values, chr, chrLength, range, sum))
        }
    }
}