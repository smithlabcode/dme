#
# Copyright (C) 2008 Cold Spring Harbor Laboratory and Andrew D Smith
# Author: Andrew D Smith
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
# USA
# 
Copyright (C) 2008 Cold Spring Harbor Laboratory and Andrew D Smith
Author: Andrew D Smith

########################################################################
========================================================================
########################################################################

This README file corresponds to the beta distribution of DME2, an
algorithm for discovering motifs de novo that best discriminate two
sets of sequences.

At the time of this writing, DME2 program works very well, and I use
it regularly. Also at the time of this writing, the DME2 algorithm is
not finished, with certain features yet to be included in the public
release. DME2 does roughly the same thing as the original DME, just
faster and more accurately. More specifically, certain parameters have
been optimized based on experience with the original DME, and both (i)
the way in which the search space is traversed, and (ii) the sets of
column types used have been modified.

DME2 has three major modes of operation: ZOOPS mode, TCM mode and a
hybrid that uses TCM for part of the search but finishes with
ZOOPS. Consult the web for general information on the meanings of
"ZOOPS" and "TCM" in the context of motif discovery. Essentially, the
ZOOPS model in DME2 will try to discover motifs that occur in (a) in a
greater number of foreground sequences than background sequences, and
(b) with a stronger top-scoring occurrence in foreground sequences
than in background sequences. So ZOOPS is "sequence-centric". TCM
tries to discover motifs with more and stronger occurrences in the
foreground than in the background, regardless of how many occur in
each sequence.

########################################################################
========================================================================
########################################################################

INSTALLING DME2

This release includes some files from the CREAD library (also
developed by me). So that CREAD does not require installation in order
to install and run DME2, the required files have been included in a
subdirectory of this archive (i.e. where this file should be located).
Hopefully the Makefile will work to build them as needed, but I
haven't tested many platforms. To compile DME2, just do the following:

$ tar -zxvf dme2_beta_2008_05_26.tgz
$ cd dme2_beta_2008_05_26
$ make

That should be all you have to do. A binary named "dme2" should be
created, and can be moved around. NOTE: currently DME2 requires the
GNU popt library to be installed and for the headers to be visible
(i.e. in some sensible location) when you try to compile
DME2. Hopefully I will replace popt sometime soon, but I like popt...

Also, although all the C++ used in DME2 is standard (AFAIK), the GCC
compilers older than ~3.4 will not compile it because they are not
compliant.

########################################################################
========================================================================
########################################################################

SELECTING BACKGROUND SEQUENCE SETS

This is the most common question I'm asked: what should I use as a
background? I used to answer that it depends on the hypothesis that
motivated the motif discovery in the first place. Nobody really likes
that answer, and I think it's because motif discovery is rarely though
of in the context of testing a hypothesis. So here are some new
answer(s):

(1) If your foreground is promoters for genes showing differential
expression between two sets of experiments, then your background might
consist of (a) promoters selected uniformly at random from the set of
all promoters, (b) promoters of genes showing the least change between
the two conditions, (c) promoters with similar single- and
di-nucleotide composition

If the same motifs come up when you use all 3 backgrounds, then you
may have found something interesting.  However, make sure you
understand why each of these backgrounds might lead to motifs that are
simply artifacts.

(2) If your foreground is regions identified with ChIP-seq (or
ChIP-chip...), then your most conservative background would be the set
of flanking regions for your ChIP-positive regions. If the
peak-calling method is not accurate, your background will be
contaminated, though. If your (ChIP-)chip only interrogates promoters,
then be very careful when using flanking sequences, as they might
differ significantly in composition from your foreground.

(3) Never simply take sequences sampled uniformly at random from the
genome!

(4) Although I always use real sequences, I often use randomly
generated sequences along with them. But make sure the random
sequences have the same higher-order nucleotide composition as the
foreground (use, e.g., the Shufflet program).

It is a good idea to use several different backgrounds.

########################################################################
========================================================================
########################################################################

ARE MY MOTIFS REAL?

I don't know. But before even considering that question, try to see if
they have significant scores. Unfortunately there is no quick and easy
way to check statistical significance of motif scores, and the reason
is probably that we don't understand real genomic (regulatory)
sequences well enough yet. So I randomly shuffle the _labels_ on the
foreground and background sequences. In other words, pool all
sequences, and select randomly some number equal to the size of the
foreground, use that as the random foreground, with the remaining as
the random background. Then run DME2 on the sequences the exact same
way as you did for your original foreground/background. Do this 100
times, and keep the scores. If your original scores are better than
any from these 100 random shuffles, you can feel good about the motifs
(but then try 10,000 shuffles before publishing). However, if there is
some bias in your background selection, this method will not help
(e.g. if G/C content differs between foreground and background).

########################################################################
========================================================================
########################################################################

RUNNING DME2

Here is how I usually run DME2:

dme2 -v -n 200 -w 10 -o output.mat -b background.fa foreground.fa

This will produce the top 200 motifs of width 10 using the "hybrid"
model. 200 might seem like lots, but currently in order for DME2 to
give good results, it requires that you ask for many motifs, even if
you are only interested in the top 3-5. A description of other
command-line options is given below. Previously, DME could only get
motifs of width up to ~12. With DME2 you can probably go up to 14 or
higher. But make sure to check the memory usage if the "-z" flag is
used.

And DME2 requires that the sequences be given in FASTA format.

########################################################################
========================================================================
########################################################################

COMMAND LINE ARGUMENTS

Usage: dme2 [OPTIONS] filename
  -z, --zoops                 use the ZOOPS model (default: hybrid)

	This flag makes the algorithm perform similar to DME-B.
	Because the run time, and therefore number of iterations of
	refinement, is so different (i.e. better), the results might
	seem qualitatively different from DME-B. With very large data
	sets memory can become a problem using this option.

  -t, --tcm                   use the TCM model (default: hybrid)

	This flag makes the algorithm more similar to the original DME,
	in terms of the objective function. This is faster than using the
	"-z" flag, but also runs the risk of finding unmasked repeats
	or motifs that occurr very frequently in a small subset of
	sequences in the foreground.

	If neither "-t" nor "-z" are used, then the identification of seeds
	is done just as with "-t", but the refinement is done as with
	"-z" (called the "hybrid" behavior). This is often the best
	balance of speed and sensitivity, which is why this behavior
	is the default.

  -b, --background=STRING     background sequence file (FASTA format)
  -o, --output=STRING         output file name (default: stdout)
  -n, --number=INT            number of motifs to produce. (default: 1)

	You should ask for many more than 1. This is different from DME,
	and this number both determines the number of seeds used in the
	search, and the number of outputs. So ask for more outputs to get
	more accuracy, because the top ranking motif when requesting 1
	motif will probably not be the same as (or as good as) the top
	when requesting 50 motifs.

  -p, --prefix=STRING         string prepended motif accession
  -w, --width=INT             minimum desired motif width (default: 8)
  -i, --bits=FLOAT            min bits per column (default depends on width)

	The value of this parameter does not affect the run time as much
	as it did in DME or DME-B.

  -c, --correction=FLOAT      correction for 0 in matrices (default: 1e-10)

	It does not make much sense to mess with this parameter too much,
	since the best performance is when this value is so low that it
	behaves like a 0, in all respects except that taking the log of it
	makes sense.

  -r, --refine=FLOAT          refinement granularity (default depends on width)

	Unlike DME, this parameter does not have much influence in DME2. Because the
	algorithm is different, it is much less sensitive to different values of this
	parameter. The default is 0.25, but 0.125 could also be tested.

  -a, --adjust=FLOAT          adjust contribution of fg and bg

	This argument is used to adjust the importance of more/stronger occurrences
	in the foreground, relative to the importance of fewer/weaker
	occurrences in the background. For example, if you think your foreground has
	much noise (i.e. many sequences that don't belong), you might want to
	_increase_ the value of this parameter so that each occurrence of the
	motif in the foreground has greater value (and hence fewer are needed). Doing
	so will generally have the effect of driving down the number of occurrences
	identified in the background for motif deemed the best by
	DME2.

  -C, --changes=INT           number of changes per refinement iteration

	This parameter defaults to 1, and it controls how "Greedy" is the
	refinement procedure. I have not used any other value for a long time
	and I find that 1 works well. Other values will take longer.

  -I, --iterations=INT        number of refinement iterations

	The default now is to iterate until convergence, which seems to work
	quite well.

  -v, --verbose               print more run information

########################################################################
========================================================================
########################################################################

Andrew D Smith
(Mon May 26 18:31:05 EDT 2008 @ CSHL)
