use strict;
use warnings;

## copied from user dongliangleng at UM

use File::Basename;

die "perl $0 <sh.list> <total_job_num>\n" unless @ARGV==2;

my ($list,$num)=@ARGV;

my $path=dirname($list);		
my $name=basename($0,".pl");
my $log="$path/$name.log";
if( ! -e $log){
        `touch $log`;
}

my $R_num=`squeue -hu 'lakhan'|grep 'NORMAL'|wc -l`;			#20

my $total_num=`cat $list|wc -l`;		#300
my $done_num=`cat $log|wc -l`;			#384
my $done_num=$done_num/3;				#384/3=128
my $left=$total_num - $done_num;		#300-128=172

my $ava=$num -$R_num;					#20-20

#print("$log\t$total_num\t$done_num\t$left\t$ava\n");
for my $a (`cat $list|tail -n $left|head -n $ava`){
        my $jobws = dirname($a);
        chdir $jobws;
        #my $cur=`pwd`;
        #print("$cur\t");
        `date '+%Y-%m-%d\ %H:%M:%S' >>$log`;
        `echo '$a' >>$log`;
        `sbatch $a`;
        chdir $path;
        #$cur=`pwd`;
        #print("$cur\n");

}

