/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package PedigreeSim;

import java.text.DecimalFormat;
import java.util.Arrays;

/**
 * Locus stores the general attributes of alocus: name, position
 * and the allelenames corresponding to each founder allele.
 * @author Roeland Voorrips
 */
public class Locus {
    private Chromosome chrom;  //on which chromosome is this locus located
    private String locusName;
    protected double position;   //in Morgan!
    private String[] alleleName; //one for each founder allele; duplicates are allowed
    private String[] uniqueAlleleNames; //only the unique allele names, in alphabetical order
    private int alleleCount=-1; //the number of unique allele names; -1 indicates: not determined yet

    public Locus(String locusName, double position, String[] alleleNames,
            PopulationData popdata) {
        if (locusName==null || locusName.trim().isEmpty()) {
            DecimalFormat fix3 = new DecimalFormat("#0.000");
            locusName = fix3.format(position);
        } else {
            this.locusName = locusName;
        }
        this.position = position;
        this.chrom = null;
        alleleName = new String[popdata.founderAlleleCount];
        if (alleleNames==null) {
            for (int i=0; i<popdata.founderAlleleCount; i++) {
                alleleName[i] = ""+i;
            }
        } else { //copy the correct number of alleleNames
            int minn = Math.min(alleleNames.length, this.alleleName.length);
            System.arraycopy(alleleNames, 0, this.alleleName, 0, minn);
            //add extra alleleNames if list not long enough:
            for (int i=minn; i<this.alleleName.length; i++) {
                this.alleleName[i] = ""+i;
            }
        }
    }

    public Chromosome getChrom() {
        return chrom;
    }

    public String getLocusName() {
        return locusName;
    }

    public double getPosition() {
        return position;
    }

    public void setChrom(Chromosome chrom) {
        this.chrom = chrom;
    }

    public String[] getAlleleNames() {
        return alleleName;
    }

    public String getAlleleName(int founder) {
        if (founder<0 || founder >=alleleName.length) {
            return "";
        } else {
            return alleleName[founder];
        }
    }

    public void setAlleleName(String[] alleleName) {
        this.alleleName = alleleName;
    }
    
    /**
     * countAlleleNames :
     * internal function, called only once, that determines the number
     * of different alleles at this locus and also the lowest and highest
     * alleleName (alphabetically). Stored in the variables
     * alleleCount, minAlleleName and maxAlleleName
     */
    private void countAlleleNames() {
        alleleCount=0;
        String[] tmp = new String[alleleName.length];
        for (int i=0; i<alleleName.length; i++) {
            int j = i-1;
            while (j>=0 && !alleleName[i].equals(alleleName[j])) {
                j--;
            }
            if (j<0) {
                //new name found
                tmp[alleleCount++] = alleleName[i];
            }
        } //for i
        uniqueAlleleNames = Arrays.copyOfRange(tmp, 0, alleleCount);
        Arrays.sort(uniqueAlleleNames);
    }

    public String getMinAlleleName() {
        if (alleleCount==-1) {
            countAlleleNames();
        }
        return uniqueAlleleNames[0];
    }

    public String getMaxAlleleName() {
        if (alleleCount==-1) {
            countAlleleNames();
        }
        return uniqueAlleleNames[alleleCount-1];
    }

    public int getAlleleCount() {
        if (alleleCount==-1) {
            countAlleleNames();
        }
        return alleleCount;
    }
    
    public String[] getUniqueAlleleNames() {
        if (alleleCount==-1) {
            countAlleleNames();
        }
        return uniqueAlleleNames;
    }

}
