/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package PedigreeSim;

import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;

/**
 * Chromosome describes the general features of a chromosome, including its
 * start- and end positions (and therefore its length) and centromere position.
 * It can also have a linkage map, implemented as an ArrayList of Locus elements.
 * The composition (in terms of founder alleles) of a specific instance
 * of a chromosome is described by a HaploStruct. The Locus alleles present
 * on such an instance follow from the Locus positions and its founder alleles.
 * @author Roeland Voorrips
 */
public class Chromosome {
    PopulationData popdata;
    String chromName;
    private double startPos;      //in Morgan; a chromosome doesn't have to start at position 0
    private double endPos;        //in Morgan
    private double centromerePos; //in Morgan
    ArrayList<Locus> locus;
    
    //Arrays for use in test mode:
    //arrays for storing the cumulative counts of each founder allele:
    //1st index=popdata.ploidy=#founder alleles, 2nd index=#markers on chrom
    int[][] founderAlleleCountAll; //over all gametes
    int[][] founderAlleleCountSel; //one random gamete per meiosis
    //array for storing the locus genotypes over all iterations
    //1st index=#locGenotypes, 2nd index=#markers on chromosome
    int[][] locGenotypeCount; //over all iterations
    //arrays for storing the recombination count between any two markers:
    //both indices 0..makercount-1
    int[][] recombCountAll; //counts recombinations over all gametes, per iteration
    int[][] recombCountSel; //counts recombinations of one random gamete per meiosis, per iteration
    //int[][] recombCountCum; //cumulates recombCountSel over all iterations
    double[][] recombFrCum; //cumulates recombination fraction over all iterations
    double[][] recombFrSqCum; //cumulates the squared recomb.fraction over all iterations
    int[][] mapdistCount; //counts the non-NaN distances
    double[][] mapdistCum; //sum of map distances over all iterations
    double[][] mapdistSqCum; //sum of squared map distances over all iterations
    //array for storing a frequency distribution of the number of
    //recombination points
    //index i is for i recombination points, length of array will be sized as necessary
    int[] recombPoints; //cumulative over all iterations
    //array for storing the frequency distribution of the number of different
    //founder alleles in one chromosome (min is 1, max is ploidy) 
    //index is count-1: 0..ploidy-1
    int[] founderCount; //cumulative over all iterations


    public Chromosome(String chromName, double startPos, double endPos,
            double centromerePos, PopulationData popdata) {
        this.popdata = popdata;
        this.chromName = chromName;
        this.centromerePos = centromerePos;
        if (centromerePos<startPos) {
            this.startPos = centromerePos; 
        } else {
            this.startPos = startPos;
        }
        if (centromerePos>endPos) {
            this.endPos = centromerePos;
        } else {
            this.endPos = endPos;
        }
        locus = new ArrayList<Locus>();
    }

    public Chromosome(String chromName, double length, double centromerePos,
            PopulationData popdata) {
        this(chromName, 0.0, length, centromerePos, popdata);
        //startPos=0.0, endPos=length
    }

    public double getCentromerePos() {
        return centromerePos;
    }

    public String getChromName() {
        return chromName;
    }

    public double getLength() {
        return endPos-startPos;
    }

    public ArrayList<Locus> getLocus() {
        return locus;
    }

    public double getStartPos() {
        return startPos;
    }

    public double getEndPos() {
        return endPos;
    }

    public PopulationData getPopdata() {
        return popdata;
    }

    public int getChromNumber() {
        return popdata.getChromosome().indexOf(this);
    }

    public void addLocus(Locus newLocus) throws Exception {
        double tolerance = 1.0e-7;
        newLocus.setChrom(this);
        if (newLocus.getPosition()<startPos) {
            if (newLocus.getPosition()<startPos-tolerance) {
                throw new Exception("addLocus: position invalid");
            }
            else {
                newLocus.position = startPos;
            }
        }
        if (newLocus.getPosition()>endPos) {
            if (newLocus.getPosition()>endPos+tolerance) {
                throw new Exception("addLocus: position invalid");
            }
            else {
                newLocus.position = endPos;
            }
        }
        if (locus.isEmpty()) locus.add(newLocus);
        else {
            int index = 0;
            Locus loc = locus.get(index);
            while (index<locus.size() && 
                    locus.get(index).getPosition()<= newLocus.getPosition()) {
                index++;
            }
            locus.add(index, newLocus);
        }
    } //addLocus

    /**
     * autoAddLoci:
     * for use in test situations where we don't want to specify a whole map;
     * not used in the regular PedigreeSim
     * @param start
     * @param end
     * @param interval
     * @param prefix 
     */
    public void autoAddLoci(double start, double end, double interval, String prefix) {
        if (prefix==null || prefix.trim().equals("")) {
            prefix = "";
        }
        if (start<startPos) start=startPos;
        if (end>endPos) end=endPos;
        int count = (int) (1.0 +0.01*interval+(end-start)/interval);
        String zeroes = "0";
        if (count>=10) zeroes += "0";
        if (count>=100) zeroes += "0";
        if (count>=1000) zeroes += "0";
        DecimalFormat df = new DecimalFormat(zeroes);
        double pos = start; count = 0;
        Locus loc; count=0;
        while (pos < end+0.01*interval) {
            loc = new Locus(prefix + df.format(count), pos, null, popdata);
            try {
                addLocus(loc);
            } catch (Exception ex) {} //cannot happen
            pos += interval;
            count++;
        }
    } //autoAddLoci

    protected int maxSegCount() {
        int maxSeg = 0;
        int chromNum = getChromNumber();
        for (Individual ind: popdata.getIndividual()) {
            for (int p=0; p<popdata.ploidy; p++) {
                int n = ind.getHaploStruct(chromNum, p).segmentCount();
                if (n>maxSeg) maxSeg=n;
            }
        }
        return maxSeg;
    } //maxSegCount

    /**
     * printMapHorizontal:
     * prints two lines, one with the marker names and one with the positions.
     * Used for generating output in test mode
     * @param out
     * @param captions 
     */
    void printMapHorizontal(PrintWriter out, boolean captions) {
        // for use in main.testMeiosis
        if (captions) out.print("marker");
        for (int loc=0; loc<getLocus().size(); loc++)
            out.print("\t"+getLocus().get(loc).getLocusName());
        out.println();
        if (captions) out.print("position");
        for (int loc=0; loc<getLocus().size(); loc++)
            out.print("\t"+getLocus().get(loc).getPosition());
        out.println();
    } //printMapHorizontal

}
