/**
 * 
 */
package ca.mcgill.genome.mps.core.parser.ace;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.lang.math.IntRange;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 *
 */
public class Read {
  private static final Logger log = LoggerFactory.getLogger(Read.class);
  private final static Pattern range = Pattern.compile("^([^\\.]+)\\.(\\d+)-(\\d+)");

  public enum Direction {
    UNKNOWN, FWD, REV
  };

  private final String name;
  private final IntRange clear;
  private int offset;
  private boolean complemented;
  private int nbPaddedBases;
  private int nbWholeReadInfoItems;
  private int nbReadTags;
  private int qualClipStart;
  private int qualClipEnd;
  private int alignClipStart;
  private int alignClipEnd;
  private Direction direction = Direction.UNKNOWN;
  private String sequence;
  private byte baseQualities[];

  public Read(String name) {
    this.name = name;
    Matcher rangeMatch = range.matcher(name);
    if (rangeMatch.find()) {
      this.clear = new IntRange(Integer.parseInt(rangeMatch.group(2)), Integer.parseInt(rangeMatch.group(3)));
    } else {
      this.clear = null;
    }
  }

  public String getRealName() {
    int idx = name.indexOf('.');
    if (idx == -1) {
      return name;
    } else {
      return name.substring(0, idx);
    }
  }

  public List<Integer> gaps() {
    List<Integer> retVal = new ArrayList<Integer>();

    String realSeq = sequence;
    if (clear != null) {
      realSeq = sequence.substring(clear.getMinimumInteger(), clear.getMaximumInteger());
    }
    for (int i = 0; i < realSeq.length() && (i + offset) > -1; i++) {
      char readBase = realSeq.charAt(i);
      if (readBase == '*') {
        retVal.add(i - retVal.size());
      }
    }
    return retVal;
  }

  public String getSequence() {
    return sequence;
  }

  public void setSequence(String sequence) {
    this.sequence = sequence;
  }

  public byte[] getBaseQualities() {
    return baseQualities;
  }

  public void setBaseQualities(byte[] baseQualities) {
    this.baseQualities = baseQualities;
  }

  public String getName() {
    return name;
  }

  public int getOffset() {
    return offset;
  }

  public void setOffset(int offset) {
    this.offset = offset;
  }

  public boolean isComplemented() {
    return complemented;
  }

  public void setComplemented(boolean complemented) {
    this.complemented = complemented;
  }

  public int getNbPaddedBases() {
    return nbPaddedBases;
  }

  public void setNbPaddedBases(int nbPaddedBases) {
    this.nbPaddedBases = nbPaddedBases;
  }

  public int getNbWholeReadInfoItems() {
    return nbWholeReadInfoItems;
  }

  public void setNbWholeReadInfoItems(int nbWholeReadInfoItems) {
    this.nbWholeReadInfoItems = nbWholeReadInfoItems;
  }

  public int getNbReadTags() {
    return nbReadTags;
  }

  public void setNbReadTags(int nbReadTags) {
    this.nbReadTags = nbReadTags;
  }

  public int getQualClipStart() {
    return qualClipStart;
  }

  public void setQualClipStart(int qualClipStart) {
    this.qualClipStart = qualClipStart;
  }

  public int getQualClipEnd() {
    return qualClipEnd;
  }

  public void setQualClipEnd(int qualClipEnd) {
    this.qualClipEnd = qualClipEnd;
  }

  public int getAlignClipStart() {
    return alignClipStart;
  }

  public void setAlignClipStart(int alignClipStart) {
    this.alignClipStart = alignClipStart;
  }

  public int getAlignClipEnd() {
    return alignClipEnd;
  }

  public void setAlignClipEnd(int alignClipEnd) {
    this.alignClipEnd = alignClipEnd;
  }

  public Direction getDirection() {
    return direction;
  }

  public void setDirection(Direction direction) {
    this.direction = direction;
  }

  public IntRange getClear() {
    return clear;
  }

}
