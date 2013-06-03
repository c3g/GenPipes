/**
 * 
 */
package ca.mcgill.genome.mps.amos.msg;

import java.io.PrintStream;
import java.math.BigDecimal;
import java.text.DecimalFormat;

/**
 * 
 */
public class Distribution implements MessageType {
  private final BigDecimal mean;
  private final BigDecimal stdDeviation;
  // DecimalFormat isn't thread safe
  private static final ThreadLocal<DecimalFormat> decimalFormatTLS = new ThreadLocal<DecimalFormat>();

  public Distribution() {
    this(BigDecimal.ZERO, BigDecimal.ZERO);
  }

  public Distribution(BigDecimal mean, BigDecimal stdDeviation) {
    this.mean = mean;
    this.stdDeviation = stdDeviation;
  }

  @Override
  public void write(PrintStream writer) {
    DecimalFormat format = getDecimalFormat();
    writer.println("{DST");
    writer.print("mea:");
    writer.println(format.format(mean));
    writer.print("std:");
    writer.println(format.format(stdDeviation));
    writer.println('}');
  }

  private static DecimalFormat getDecimalFormat() {
    DecimalFormat retVal = decimalFormatTLS.get();
    if (retVal == null) {
      retVal = new DecimalFormat("#######0.000");
      decimalFormatTLS.set(retVal);
    }
    return retVal;
  }
}
