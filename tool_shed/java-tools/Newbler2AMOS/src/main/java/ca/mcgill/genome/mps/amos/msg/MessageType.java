package ca.mcgill.genome.mps.amos.msg;

import java.io.IOException;
import java.io.PrintStream;

public interface MessageType {
  public static final int DEFAULT_LINE_LENGTH = 70;
  public static final int DEFAULT_QUALITY_SCALE = 48;

  void write(PrintStream writer) throws IOException;
}
