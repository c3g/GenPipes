/**
 * 
 */
package ca.mcgill.genome.mps.core.parser;

import org.apache.commons.lang.math.IntRange;

/**
 * 
 */
public class TrimStatus {
  private final String accno;
  private IntRange trimpoints;
  private Integer trimmedLength;
  private IntRange origTrimpoints;
  private Integer origTrimmedLength;
  private Integer rawLength;

  public TrimStatus(String accno) {
    this.accno = accno;
  }

  public String getAccno() {
    return accno;
  }

  public IntRange getTrimpoints() {
    return trimpoints;
  }

  public void setTrimpoints(IntRange trimpoints) {
    this.trimpoints = trimpoints;
  }

  public Integer getTrimmedLength() {
    return trimmedLength;
  }

  public void setTrimmedLength(Integer trimmedLength) {
    this.trimmedLength = trimmedLength;
  }

  public IntRange getOrigTrimpoints() {
    return origTrimpoints;
  }

  public void setOrigTrimpoints(IntRange origTrimpoints) {
    this.origTrimpoints = origTrimpoints;
  }

  public Integer getOrigTrimmedLength() {
    return origTrimmedLength;
  }

  public void setOrigTrimmedLength(Integer origTrimmedLength) {
    this.origTrimmedLength = origTrimmedLength;
  }

  public Integer getRawLength() {
    return rawLength;
  }

  public void setRawLength(Integer rawLength) {
    this.rawLength = rawLength;
  }
}
