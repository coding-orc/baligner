/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package core.sb;

/**
 *
 * @author jarcher
 */
public class Seq{
    public int i;
    public String aligned, voucher, accession, dna;
    public String key;
    public boolean allNs = false;
    
    public Seq(String key, String dna){
        this.key = key;
        this.dna = dna;
        allNs = false;
    }
    
    public Seq(String aligned, String voucher, String accession, String dna, int i){
        this.aligned = aligned;
        this.voucher = voucher;
        this.accession = accession;
        this.dna = dna;
        this.i = i;
        allNs = false;
    }
}