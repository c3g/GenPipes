/*
 * Copyright 2010© McGill University and Génome Québec Innovation Centre. 
 * All rights reserved.
 *
 */

package ca.mcgill.genome.mps;

import java.io.File;
import java.math.BigDecimal;

/**
 *
 */
public class SffFileDetail {
  private final File sffFile;
  private BigDecimal insertMean = null;
  private BigDecimal insertDev = null;

  public SffFileDetail(File sffFile) {
    this.sffFile = sffFile;
  }

  public BigDecimal getInsertMean() {
    return insertMean;
  }

  public void setInsertMean(BigDecimal insertMean) {
    this.insertMean = insertMean;
  }

  public BigDecimal getInsertDev() {
    return insertDev;
  }

  public void setInsertDev(BigDecimal insertDev) {
    this.insertDev = insertDev;
  }

  public File getSffFile() {
    return sffFile;
  }

}
