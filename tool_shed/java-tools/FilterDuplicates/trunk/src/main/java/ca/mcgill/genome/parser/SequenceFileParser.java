package ca.mcgill.genome.parser;

import java.io.IOException;

public abstract class SequenceFileParser {

  public abstract Sequence nextSequence() throws IOException;

  /**
   * This method doesn't guarantee to be precise.
   * @return
   */
  public abstract long getByteCount();

  public abstract void close() throws IOException;
}
