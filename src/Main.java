

import com.github.SnowFlakes.File.GffFile.GFF3File;
import com.github.SnowFlakes.IO.GFF3ReaderExtension;
import com.github.SnowFlakes.Utils.PetCluster;
import com.github.SnowFlakes.unit.ChrRegion;
import com.github.SnowFlakes.unit.Region;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.gff.Gff3Feature;
import org.apache.commons.math3.util.MathUtils;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

public class Main {
    public File bed_graph_file;
    public File annotation_file;
    public File out_file;

    public static void main(String[] args) {
//        args = new String[]{"s9am1.align.RPKM.bedgraph", "XENTR_10.0_GCF.Chr2.gff3"};
        if (args.length<3){
            System.err.println("Usage: java -jar CoverageCalculation.jar <bedgraph file> <gff3 file> <out file>");
            System.exit(1);
        }
        Main main = new Main();
        main.bed_graph_file = new File(args[0]);
        main.annotation_file = new File(args[1]);
        main.out_file=new File(args[2]);
        try {
            main.run();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void run() throws IOException {
        System.err.println("read gff file");
        GFF3ReaderExtension gff3_reader = new GFF3File(annotation_file.getPath()).getReader();
        Gff3Feature feature;
        HashMap<String, ArrayList<ExonList>> gene_list = new HashMap<>();
        while ((feature = gff3_reader.ReadRecord()) != null) {
            if (!feature.isTopLevelFeature()) {
                continue;
            }
            if (!gene_list.containsKey(feature.getContig())) {
                gene_list.put(feature.getContig(), new ArrayList<>());
            }
            ExonList gene = new ExonList(new ChrRegion(feature.getContig(), feature.getStart(), feature.getEnd(), feature.getStrand()));
            gene.name = feature.getName();
            for (Gff3Feature f : feature.getChildren()) {
                if (f.getChildren() == null && f.getType().compareToIgnoreCase("exon") == 0) {
                    gene.exon_list.add(new Region(Math.min(f.getStart(), f.getEnd()), Math.max(f.getStart(), f.getEnd())));
                } else {
                    for (Gff3Feature ff : f.getChildren()) {
                        if (ff.getType().compareToIgnoreCase("exon") == 0) {
                            gene.exon_list.add(new Region(Math.min(ff.getStart(), ff.getEnd()), Math.max(ff.getStart(), ff.getEnd())));
                        }
                    }
                }
            }
            Collections.sort(gene.exon_list);
            gene.exon_merge();
            int sum = 0;
            for (Region r : gene.exon_list) {
                sum += r.getLength() + 1;
            }
            gene.pos_coverage = new float[sum];
            gene_list.get(gene.gene.Chr).add(gene);
        }
        gff3_reader.close();
        for (String key : gene_list.keySet()) {
            Collections.sort(gene_list.get(key));
        }
        System.err.println("read bedgraph file");
        BufferedReader bg_reader = IOUtil.openFileForBufferedReading(bed_graph_file);
        String lines;
        while ((lines= bg_reader.readLine())!=null){
            String[] strs = lines.split("\\s+");
            ChrRegion pos = new ChrRegion(strs);
            float value = Float.parseFloat(strs[3]);
            if (gene_list.containsKey(pos.Chr)){
                for(ExonList l: gene_list.get(pos.Chr)){
                    if (pos.IsOverlap(l.gene)){
                        l.add_pos(pos.region, value);
                    }else if (l.gene.region.Start>pos.region.End){
                        break;
                    }
                }
            }
        }
        BufferedWriter writer = new BufferedWriter(new FileWriter(out_file));
        int interval=100;
        for (String key: gene_list.keySet()){
            for (ExonList e: gene_list.get(key)){
                if (e.pos_coverage.length<interval){
                    continue;
                }
                float len=(float) e.pos_coverage.length/interval;
                writer.write(e.name+" ");
                if (e.gene.Orientation== Strand.REVERSE){
                    for (int i = interval-1; i >=0 ; i--) {
                        writer.write(means(e.pos_coverage,(int) len*i,(int) len*(i+1))+" ");
                    }
                }else {
                    for (int i = 0; i <interval ; i++) {
                        writer.write(means(e.pos_coverage,(int) len*i,(int) len*(i+1))+" ");
                    }
                }
                writer.write(e.pos_coverage.length+"\n");
            }
        }
        writer.close();
    }
    private float means(float[] a,int i,int j){
        float sum=0;
        for (int k = i; k <j ; k++) {
            sum+=a[k];
        }
        return sum/(j-i);
    }
}

class ExonList implements Comparable<ExonList> {
    public ChrRegion gene;
    public String name;
    public ArrayList<Region> exon_list = new ArrayList<>();
    public float[] pos_coverage;

    public ExonList(ChrRegion c) {
        gene = c;
    }

    @Override
    public int compareTo(ExonList o) {
        return this.gene.compareTo(o.gene);
    }

    /**
     * 将所有外显子合并
     */
    public void exon_merge() {
        int OldSize = 0;
        while (OldSize != exon_list.size()) {
            OldSize = exon_list.size();
            for (int i = 0; i < exon_list.size(); ++i) {
                for (int j = i + 1; j < exon_list.size(); ++j) {
                    if (exon_list.get(i).End < exon_list.get(j).Start) {
                        break;
                    }
                    exon_list.get(i).End = Math.max(exon_list.get(i).End, exon_list.get(j).End);
                    exon_list.remove(j);
                    j--;
                }
            }
        }
    }

    public void add_pos(Region r,float v) {
        if (pos_coverage == null) {
            return;
        }
        int index=0;
        for (Region s : exon_list) {
            if (s.Start>r.End) {
                break;
            }
            if (s.IsOverlap(r)){
                int s_index=Math.max(0,r.Start-s.Start)+index;
                int e_index=Math.min(r.End-s.End,0)+index+s.getLength()+1;
                for (int i = s_index; i <e_index ; i++) {
                    pos_coverage[i]+=v;
                }
            }
            index+=s.getLength()+1;
        }
    }
}