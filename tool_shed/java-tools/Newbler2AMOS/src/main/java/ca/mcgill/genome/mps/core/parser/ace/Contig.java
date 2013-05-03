/**
 * 
 */
package ca.mcgill.genome.mps.core.parser.ace;

import java.io.ByteArrayOutputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 *
 */
public class Contig {
  private static final Logger log = LoggerFactory.getLogger(Contig.class);

  private final String name;
  private String unTrimmedSequence = null;
  private byte unTrimmedBaseQualities[] = null;
  private String sequence;
  private byte baseQualities[];
  private int numBases;
  private int numColumns;
  private int numReads;
  private int numBaseSegs;
  private boolean complemented;
  private final Map<String, Read> reads = new HashMap<String, Read>();
  private final List<String> bsLines = new ArrayList<String>();

  public Contig(String name) {
    this.name = name;
  }

  public void addRead(Read read) {
    Read old = reads.put(read.getName(), read);
    if (old != null) {
      log.error("Read '{}' already set on contig '{}'", read.getName(), name);
    }
  }

  public void setUnTrimmedSequence(String unTrimmedSequence) {
    this.unTrimmedSequence = unTrimmedSequence;
    computeTrimmed();
  }

  public void setUnTrimmedBaseQualities(byte[] unTrimmedBaseQualities) {
    this.unTrimmedBaseQualities = unTrimmedBaseQualities;
    computeTrimmed();
  }

  private void computeTrimmed() {
    if (sequence != null && baseQualities != null) {
      throw new IllegalStateException("Sequence already computed");
    }

    StringBuilder trimmedSeq = new StringBuilder();
    ByteArrayOutputStream trimmedBases = new ByteArrayOutputStream();
    ByteArrayOutputStream gappedUnTrimmedBaseQualities = new ByteArrayOutputStream();
    if (unTrimmedBaseQualities != null && unTrimmedSequence != null) {
      for (int seqIdx = 0, qualIdx = 0; seqIdx < unTrimmedSequence.length(); seqIdx++) {
        if (unTrimmedSequence.charAt(seqIdx) != '*') {
          trimmedSeq.append(unTrimmedSequence.charAt(seqIdx));
          trimmedBases.write(unTrimmedBaseQualities[qualIdx]);
          gappedUnTrimmedBaseQualities.write(unTrimmedBaseQualities[qualIdx]);
          qualIdx++;
        } else {
          gappedUnTrimmedBaseQualities.write((byte) 1);
        }
      }
      unTrimmedBaseQualities = gappedUnTrimmedBaseQualities.toByteArray();
      sequence = trimmedSeq.toString();
      baseQualities = trimmedBases.toByteArray();
    }
  }

  public int getNumBases() {
    return numBases;
  }

  public void setNumBases(int numBases) {
    this.numBases = numBases;
  }

  public int getNumColumns() {
    return numColumns;
  }

  public void setNumColumns(int numColumns) {
    this.numColumns = numColumns;
  }

  public int getNumReads() {
    return numReads;
  }

  public void setNumReads(int numReads) {
    this.numReads = numReads;
  }

  public int getNumBaseSegs() {
    return numBaseSegs;
  }

  public void setNumBaseSegs(int numBaseSegs) {
    this.numBaseSegs = numBaseSegs;
  }

  public String getUnTrimmedSequence() {
    return unTrimmedSequence;
  }

  public byte[] getUnTrimmedBaseQualities() {
    return unTrimmedBaseQualities;
  }

  public String getSequence() {
    return sequence;
  }

  public byte[] getBaseQualities() {
    return baseQualities;
  }

  public boolean isComplemented() {
    return complemented;
  }

  public void setComplemented(boolean complemented) {
    this.complemented = complemented;
  }

  public String getName() {
    return name;
  }

  public Map<String, Read> getReads() {
    return reads;
  }

  public List<String> getBsLines() {
    return bsLines;
  }
}
