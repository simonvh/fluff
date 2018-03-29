import pytest
import os
import urllib.request
import tarfile
from tempfile import NamedTemporaryFile
from sklearn.metrics import v_measure_score

from fluff.parse import parse_cmds
from fluff.commands.profile import profile
from fluff.commands.heatmap import heatmap
from fluff.commands.bandplot import bandplot

@pytest.fixture
def bamfile():
    return "tests/data/H3K4me3.bam"

@pytest.fixture
def bwfile():
    return "tests/data/profile.bw"

@pytest.fixture
def bedfile():
    return "tests/data/profile.bed"

@pytest.fixture
def regionfile():
    return "tests/data/profile_region.bed"

@pytest.fixture
def test_data_from_osf():
    fnames = [
        "tests/data/big/H1_H3K27ac.bam",
        "tests/data/big/H1_H3K27ac.bam.bai",
        "tests/data/big/mesenchymal_H3K27ac.bam",
        "tests/data/big/mesenchymal_H3K27ac.bam.bai",
        "tests/data/big/mesendoderm_H3K27ac.bam",
        "tests/data/big/mesendoderm_H3K27ac.bam.bai",
        "tests/data/big/neuronal_progenitor_H3K27ac.bam",
        "tests/data/big/neuronal_progenitor_H3K27ac.bam.bai",
        "tests/data/big/trophoblast_H3K27ac.bam",
        "tests/data/big/trophoblast_H3K27ac.bam.bai",
        "tests/data/big/peaks.bed",
        ]
    
    download = False
    for fname in fnames:
        if not os.path.exists(fname):
            download = True
            break

    if download:
        # test data tarball on osf.io
        url = "https://osf.io/6yftg/download"
        tarball = "tests/data/big/test_data.tgz"
        urllib.request.urlretrieve(url, tarball)
        with tarfile.open(tarball) as tf:
            tf.extractall(path="tests/data/big/")
        os.unlink(tarball)
    
    clusters = "tests/data/big/clusters.kmeans.euclidean.5.txt"
    return fnames[-1], [f for f in fnames if f[-3:] == "bam"], clusters

def test_profile(bamfile):
    # Only tests of the command runs successfully,
    # doesnt't check the image
    with NamedTemporaryFile(prefix="fluff.", suffix=".png") as tmp:
        args = ["profile",
                "-i", "scaffold_1:44749422-44750067",
                "-d", bamfile,
                "-o", tmp.name]
        args = parse_cmds().parse_args(args)
        profile(args)

def test_heatmap(bamfile, bedfile, bwfile, regionfile):
    # Only tests of the command runs successfully,
    # doesnt't check the image
    with NamedTemporaryFile(prefix="fluff.", suffix=".png") as tmp:
        args = ["heatmap",
            "-f", regionfile,
            "-d", bamfile, bwfile, bedfile,
            "-o", tmp.name]
        args = parse_cmds().parse_args(args)
        heatmap(args)

def test_plots_big(test_data_from_osf):
    peaks, bamfiles, clusters = test_data_from_osf
    with NamedTemporaryFile(prefix="fluff.", suffix=".png") as f:
        args = [
                "heatmap", 
                '-f', peaks,
                "-d", *bamfiles,
                "-o", f.name,
                "-C", "kmeans",
                "-k", "5",
                ]
        args = parse_cmds().parse_args(args)
        heatmap(args)
        
        # Reading current clusters
        fcluster = f.name + "_clusters.bed"
        pred_clusters = []
        for line in open(fcluster):
            vals = line.strip().split("\t")
            pred_clusters.append([f"{vals[0]}:{vals[1]}-{vals[2]}", vals[4]])
        
        # Reading reference clusters
        cmp_clusters = []
        for line in open(clusters):
            vals = line.strip().split("\t")
            cmp_clusters.append(vals)

        # sort by region name
        cmp_clusters = [x[1] for x in sorted(cmp_clusters)]
        pred_clusters = [x[1] for x in sorted(pred_clusters)]
        
        v = v_measure_score(cmp_clusters, pred_clusters)

        assert v > 0.25

        # Test bandplot
        args = [
                "bandplot", 
                '-f', fcluster,
                "-d", *bamfiles,
                "-o", f.name,
                ]
        args = parse_cmds().parse_args(args)
        bandplot(args)
    
        tmpnames = [
                f.name + "_clusters.bed",
                f.name + "_readCounts.txt",
                ]

        for name in tmpnames:
            os.unlink(name)
