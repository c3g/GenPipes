/*
 * Copyright 2010© McGill University and Génome Québec Innovation Centre. 
 * All rights reserved.
 *
 */

package ca.mcgill.genome.mps.amos.msg;

import java.io.IOException;
import java.io.PrintStream;
import java.util.List;

/**
 *
 */
public class Tile implements MessageType {
  private final long src;
  private final int clearMin;
  private final int clearMax;
  private final int offset;
  private final List<Integer> gaps;

  public Tile(long src, int clearMin, int clearMax, int offset, List<Integer> gaps) {
    this.src = src;
    this.clearMin = clearMin;
    this.clearMax = clearMax;
    this.offset = offset;
    this.gaps = gaps;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void write(PrintStream writer) throws IOException {
    writer.println("{TLE");
    writer.print("src:");
    writer.println(src);
    writer.print("off:");
    writer.println(offset);
    writer.print("clr:");
    writer.print(clearMin);
    writer.print(',');
    writer.println(clearMax);
    if (gaps != null) {
      writer.println("gap:");
      boolean first = true;
      for (Integer gap : gaps) {
        if (!first) {
          writer.print(' ');
        } else {
          first = false;
        }
        writer.print(gap);
      }
      writer.println();
      writer.println('.');
    }
    writer.println('}');
  }

  public long getSrc() {
    return src;
  }

  public int getClearMin() {
    return clearMin;
  }

  public int getClearMax() {
    return clearMax;
  }

  public int getOffset() {
    return offset;
  }

  public List<Integer> getGaps() {
    return gaps;
  }

}
