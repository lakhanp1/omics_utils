
```bash
echo $0
#-bash
```

```bash
sh
#sh-3.2$
#sh-3.2$ echo $0
#sh
#sh-3.2$ bash
```

```bash
echo $0
#bash
```

```bash
echo $SHELL
#/bin/bash
```

```bash
ls /etc/host?
#/etc/hosts
```

```bash
ls /etc/host?.???
#ls: /etc/host?.???: No such file or directory
```

```bash
ls /etc/host?.????
#/etc/hosts.deny
```

```bash
ls /etc/host?.[ad]
#ls: /etc/host?.[ad]: No such file or directory
```

```bash
ls /etc/host?.[ad]*
#/etc/hosts.allow  /etc/hosts.deny
```

```bash
ls /etc/[r-w]*.conf
#/etc/radvd.conf   /etc/resolv.conf        /etc/sensors.conf   /etc/smartd.conf               /etc/sysctl.conf  /etc/updatedb.conf   /etc/webalizer.conf
#/etc/reader.conf  /etc/scrollkeeper.conf  /etc/sestatus.conf  /etc/smt_dhcp_ra_support.conf  /etc/syslog.conf  /etc/warnquota.conf  /etc/wvdial.conf
```

```bash
date
#Mon Jun 15 10:32:54 IST 2015
```

```bash
#joining more than one command
```

```bash
date; whoami; echo hello
#Mon Jun 15 10:33:45 IST 2015
#bashshell14
#hello
```

#### making one command depend on other command

```bash
cd /etc && ls -l passwd
#-rw-r--r-- 1 root root 48732 Jun 12 16:41 passwd
```

```bash
cd
cd /etccc && ls -l passwd
#bash: cd: /etccc: No such file or directory
```

#### correcting the previous command

```bash
^etccc^etc
#bash: :s^etccc^etc: substitution failed
```

```bash
touch /etc/foo
#touch: cannot touch `/etc/foo`: Permission denied
```

```bash
touch /etc/foo || echo you are not root
#touch: cannot touch `/etc/foo`: Permission denied
#you are not root
```

## shell variables

```bash
a=foo
```

```bash
echo $a
#foo
```

#### rule1: when assigning a value to the var, NEVER give space before and after =

```bash
a=hello
```

```bash
echo $a
#hello
```

```bash
a=date
#what will this produce?
```

```bash
echo $a
#date

#we did not want this, did we?

```

```bash
a=`date`
```

```bash
echo $a
#Mon Jun 15 10:39:44 IST 2015

#the command is evaluated 1st and the result is assigned to var
#this is called as COMMAND SUBSTITUTION
```

```bash
#can be also written as,
a=$(pwd)
```

```bash
echo $a
#/home/bashshell14
```

```bash
which chmod
#/bin/chmod
```

```bash
ls -l
#total 0
```

```bash
ls -l | which chmod
#/bin/chmod
```

## command redirection

```bash
ls -l `which chmod`
#-rwxr-xr-x 1 root root 38564 Feb 23  2010 /bin/chmod
```

```bash
a=foo
echo $a
#foo
```

```bash
(a=bar)
echo $a
#foo
```

```bash
(a=bar; echo $a)
#bar
```

```bash
(a=bar; echo $a); echo $a
#bar
#foo
```

```bash
a=foo
echo $a
#foo
```

```bash
echo ${a}l
#fool
```

```bash
echo hello world
#hello world
```

```bash
echo hello      world
#hello world

#this is because shell ignores the additional spaces

echo 'hello     world'
#hello     world
```

```bash
echo "hello     world"
#hello     world

echo hello \ \ \ \ \ world
#hello      world
```

```bash
a=foo

echo "$a"
#foo

echo '$a'
#$a

echo \$a
#$a
```

#### rule2: whenever in doubt, use '' except when its $

```bash
a='hello     world'
echo $a
#hello world

echo "$a"
#hello     world
```

## shell arithematic

```bash
a=5
b=7
c=$a+$b
echo $c
#5+7
```

```bash
c=`$a+$b`
#bash: 5+7: command not found
```

```bash
expr 6 + 9
#15

c=`expr $a + $b`

echo $c
#12
```

```bash
d=`expr $a * $b`

echo $\d
#$d

echo $d
#35
```

```bash
ls

#since there are no files in the current dir, the * is working as multiplication
#if a file is present then * will extrapolate to the files and it will fail

vi tmp.txt

expr 5 * 7
#expr: syntax error

#this is not working because * is extrapolating to file name tmp.txt
#so need to escape it

expr 5 \* 7
#35
```

```bash
echo my lucky *
#my lucky tmp.txt

echo my lucky \*
#my lucky *
```

```bash
touch H

echo whats up ?
#whats up H

echo 'whats up ?'
#whats up ?

echo whats up \?
#whats up ?
```

```bash
let d=a*b

echo $d
#35

((e=a*b))
echo $e
#35
```

```bash
#let command tells shell that this is an arithematic operation
let d=a*c
echo $d
#60
```

## arrays in shell

```bash
arr=(foo bar car)
```

```bash
arr=(foo bar car persistent)
```

```bash
echo ${arr[*]}
#foo bar car persistent
```

```bash
echo ${#arr[*]}
#4
```

```bash
echo ${arr[3]}
#persistent
```

```bash
echo ${#arr[3]}
#10
```

```bash
arr+=(amdocs ibm)
```

```bash
echo ${arr[*]}
#foo bar car persistent amdocs ibm
```

```bash
fname=(/etc/*.conf)
```

```bash
echo ${fname[*]}
#/etc/ant.conf /etc/autofs_ldap_auth.conf /etc/capi.conf /etc/cdrecord.conf /etc/conman.conf /etc/cyrus.conf /etc/dhcp6c.conf /etc/dhcp6s.conf /etc/dhcpd.conf /etc/dnsmasq.conf /etc/esd.conf /etc/gpm-root.conf /etc/grub.conf /etc/gssapi_mech.conf /etc/hba.conf /etc/host.conf /etc/idmapd.conf /etc/imapd.conf /etc/initlog.conf /etc/jwhois.conf /etc/kdump.conf /etc/krb5.conf /etc/ldap.conf /etc/ld.so.conf /etc/lftp.conf /etc/libaudit.conf /etc/libuser.conf /etc/logrotate.conf /etc/ltrace.conf /etc/marimba.conf /etc/mke2fs.conf /etc/modprobe.conf /etc/mtools.conf /etc/multipath.conf /etc/nscd.conf /etc/nsswitch.conf /etc/ntp.conf /etc/pam_smb.conf /etc/pear.conf /etc/prelink.conf /etc/radvd.conf /etc/reader.conf /etc/resolv.conf /etc/scrollkeeper.conf /etc/sensors.conf /etc/sestatus.conf /etc/smartd.conf /etc/smt_dhcp_ra_support.conf /etc/sysctl.conf /etc/syslog.conf /etc/updatedb.conf /etc/warnquota.conf /etc/webalizer.conf /etc/wvdial.conf /etc/xinetd.conf /etc/yp.conf /etc/ypserv.conf /etc/yum.conf
```

```bash
unset arr[1]
```

```bash
echo ${arr[*]}
#foo car persistent amdocs ibm
```

```bash
echo "${arr[*]}"
#foo car persistent amdocs ibm
```

```bash
echo '${arr[*]}'
#${arr[*]}
```

```bash
unset arr[1]
```

```bash
echo ${arr[*]}
#foo car persistent amdocs ibm
```

```bash
echo "${arr[*]}"
#foo car persistent amdocs ibm
```

```bash
echo '${arr[*]}'
#${arr[*]}
```

```bash
echo ${#arr[*]}
#5
```

```bash
a=foo
```

```bash
echo $a
#foo
```

```bash
unset a
```

```bash
echo $a
#
```

```bash
read a
#foo
```

```bash
echo $a
#foo
```

```bash
read a b
#foo bar
```

```bash
echo $a
#foo
```

```bash
echo $b
#bar
```

```bash
read a b c
#foo bar car tar
```

```bash
echo $a
#foo
```

```bash
echo $b
#bar
```

```bash
echo $c
#car tar
```

```bash
read a b
#"foo bar" car tar
```

```bash
echo $a
#"foo
```

```bash
echo $b
#bar" car tar
```

```bash
#shortcuts
```

```bash
echo the seven rules
#the seven rules
```

```bash
echo rainbow has !!:2 colors
#echo rainbow has seven colors
#rainbow has seven colors
```

```bash
!!
#echo rainbow has seven colors
#rainbow has seven colors
```

```bash
!da
#date; whoami; echo hello
#Mon Jun 15 12:55:17 IST 2015
#bashshell14
#hello
```

## Aliases for commands

```bash
alias c=clear
```

```bash
c
```

```bash
alias d='date -u'
```

```bash
d
#Mon Jun 15 07:27:51 UTC 2015
```

```bash
alias l='ls -a -l'
```

```bash
l
#total 64
#drwx------   4 bashshell14 bashshell14  4096 Jun 15 11:12 .
#drwxrwxrwx 616 root        root        20480 Jun 12 16:41 ..
#-rw-r--r--   1 bashshell14 bashshell14    33 Jun 12 16:40 .bash_logout
#-rw-r--r--   1 bashshell14 bashshell14   176 Jun 12 16:40 .bash_profile
#-rw-r--r--   1 bashshell14 bashshell14   124 Jun 12 16:40 .bashrc
#-rw-r--r--   1 bashshell14 bashshell14   515 Jun 12 16:40 .emacs
#-rw-rw-r--   1 bashshell14 bashshell14     0 Jun 15 11:12 H
#drwxr-xr-x   3 bashshell14 bashshell14  4096 Jun 12 16:40 .kde
#drwxr-xr-x   4 bashshell14 bashshell14  4096 Jun 12 16:40 .mozilla
#-rw-rw-r--   1 bashshell14 bashshell14     5 Jun 15 11:09 tmp.txt
#-rw-------   1 bashshell14 bashshell14   593 Jun 15 11:09 .viminfo
#-rw-r--r--   1 bashshell14 bashshell14   658 Jun 12 16:40 .zshrc
```

```bash
alias
#alias c='clear'
#alias d='date -u'
#alias l='ls -a -l'
#alias l.='ls -d .* --color=tty'
#alias ll='ls -l --color=tty'
#alias ls='ls --color=tty'
#alias vi='vim'
#alias which='alias | /usr/bin/which --tty-only --read-alias --show-dot --show-tilde'
```

```bash
alias d='date'
```

```bash
d -u
#Mon Jun 15 07:29:54 UTC 2015
```

```bash
unalias d
```

```bash
unalias c
```

```bash
unalias l
```

```bash
alias
#alias l.='ls -d .* --color=tty'
#alias ll='ls -l --color=tty'
#alias ls='ls --color=tty'
#alias vi='vim'
#alias which='alias | /usr/bin/which --tty-only --read-alias --show-dot --show-tilde'
```

```bash
alias pwd=date
```

```bash
pwd
#Mon Jun 15 13:01:47 IST 2015
```

```bash
\pwd
#/home/bashshell14
```

```bash
ifconfig
#bash: ifconfig: command not found
```

```bash
/sbin/ifconfig
#eth1      Link encap:Ethernet  HWaddr 00:1B:11:15:E5:D8
#          inet addr:10.44.204.213  Bcast:10.44.207.255  Mask:255.255.252.0
#          inet6 addr: fe80::21b:11ff:fe15:e5d8/64 Scope:Link
#          UP BROADCAST RUNNING MULTICAST  MTU:1500  Metric:1
#          RX packets:3758871 errors:0 dropped:0 overruns:0 frame:0
#          TX packets:211651 errors:0 dropped:0 overruns:0 carrier:0
#          collisions:0 txqueuelen:1000
#          RX bytes:290284859 (276.8 MiB)  TX bytes:22683777 (21.6 MiB)
#          Interrupt:217 Base address:0xef00
#
#lo        Link encap:Local Loopback
#          inet addr:127.0.0.1  Mask:255.0.0.0
#          inet6 addr: ::1/128 Scope:Host
#          UP LOOPBACK RUNNING  MTU:16436  Metric:1
#          RX packets:418356 errors:0 dropped:0 overruns:0 frame:0
#          TX packets:418356 errors:0 dropped:0 overruns:0 carrier:0
#          collisions:0 txqueuelen:0
#          RX bytes:46826353 (44.6 MiB)  TX bytes:46826353 (44.6 MiB)
#
```

```bash
alias ifconfig=/sbin/ifconfig
```

```bash
ifconfig
#eth1      Link encap:Ethernet  HWaddr 00:1B:11:15:E5:D8
#          inet addr:10.44.204.213  Bcast:10.44.207.255  Mask:255.255.252.0
#          inet6 addr: fe80::21b:11ff:fe15:e5d8/64 Scope:Link
#          UP BROADCAST RUNNING MULTICAST  MTU:1500  Metric:1
#          RX packets:3759996 errors:0 dropped:0 overruns:0 frame:0
#          TX packets:212362 errors:0 dropped:0 overruns:0 carrier:0
#          collisions:0 txqueuelen:1000
#          RX bytes:290384664 (276.9 MiB)  TX bytes:22769763 (21.7 MiB)
#          Interrupt:217 Base address:0xef00
#
#lo        Link encap:Local Loopback
#          inet addr:127.0.0.1  Mask:255.0.0.0
#          inet6 addr: ::1/128 Scope:Host
#          UP LOOPBACK RUNNING  MTU:16436  Metric:1
#          RX packets:418356 errors:0 dropped:0 overruns:0 frame:0
#          TX packets:418356 errors:0 dropped:0 overruns:0 carrier:0
#          collisions:0 txqueuelen:0
#          RX bytes:46826353 (44.6 MiB)  TX bytes:46826353 (44.6 MiB)
#
```

```bash
foo(){
#> echo I am from foo
#> echo hi world
#> date
#> whoami
#> }
```

```bash
foo
#I am from foo
#hi world
#Mon Jun 15 13:28:28 IST 2015
#bashshell14
```

```bash
typeset -f
#foo ()
#{
#    echo I am from foo;
#    echo hi world;
#    date;
#    whoami
#}
```

```bash
typeset -F
#declare -f foo
```

```bash
unset -f foo
```

```bash
#order to check the variable: aliases -> functions -> shell builtings -> PATH
```

```bash
vim fun_file
```

```bash
cat fun_file
#foo(){
#        echo I am from foo
#}
#
#foo2(){
#        echo U are from foo2
#}
```

```bash
source fun_file
```

```bash
foo
#I am from foo
```

```bash
foo2
#U are from foo2
```

```bash
. ./fun_file
```

```bash
#make commands/aliases/vars permanant: make changes in /etc/profile file which is owned by root
```

```bash
#since this file is owned by root we cannot make changes to it
```

```bash
#so use .bash_profile file
```

```bash
vim .bash_profile
```

```bash
#su command
```

```bash
whoami
#bashshell14
```

```bash
su bashshell17
#Password:
```

```bash
whoami
#bashshell17
```

```bash
exit
#exit
```

```bash
whoami
#bashshell14
```

```bash
a=hello
```

```bash
echo ${#a}
#5
```

```bash
echo ${a#h}
#ello
```

```bash
echo ${a#??}
#llo
```

```bash
echo ${a%??}
#hel
```

```bash
#show value of var b if its set else show the string provided
```

```bash
echo ${g-lakhan}
#lakhan
```

```bash
echo ${a:3}
#lo
```

```bash
echo ${a:1:3}
#ell
```

```bash
date
#Mon Jun 15 15:49:00 IST 2015
```

```bash
echo $?
#0
```

```bash
dateee
#bash: dateee: command not found
```

```bash
echo $?
#127
```

```bash
date ---y
#date: unrecognized option `---y`
#Try `date --help` for more information.
```

```bash
echo $?
#1
```

```bash
#redirection
```

```bash
cat
#hello
#hello
#world
#world
```

```bash
wc -l
#hello
#world
#2
```

```bash
cat nnn
#cat: nnn: No such file or directory
```

```bash
cat < nnn
#-bash: nnn: No such file or directory
```

```bash
#in first case cat gave error and in 2nd case bash gave error. this is bcoz shell is opening the file for reading
```

```bash
wc -l /etc/host
#wc: /etc/host: No such file or directory
```

```bash
wc -l /etc/hosts
#5 /etc/hosts
```

```bash
wc -l < /etc/hosts
#5
```

```bash
#Rule3: in case of <,> or >>, the RHS is always a file and in case of >, the file is truncated to 0 bytes
```

```bash
date > foo
```

```bash
date > foo
```

```bash
date > foo
```

```bash
cat foo
#Mon Jun 15 16:23:35 IST 2015
```

```bash
set -o noclobber
```

```bash
#this wont allow to truncate the file
```

```bash
date > foo
#-bash: foo: cannot overwrite existing file
```

```bash
#this is why we get the error as > is trying to truncate the file
```

```bash
#to reset it
```

```bash
set +o noclobber
```

```bash
date > foo
```

```bash
#now working
```

```bash
#to just truncate the file
```

```bash
> foo
```

```bash
: > foo # alternative way
```

```bash
cat foo
```

```bash
echo ${a}
#
```

```bash
echo ${a-word}
#word
```

```bash
a=
```

```bash
echo ${a-word}
#
```

```bash
echo ${a:-word}
#word
```

```bash
date > foo > bar
```

```bash
cat foo
```

```bash
cat bar
#Tue Jun 16 09:48:56 IST 2015
```

```bash
# nothing is written in foo file. the output of date command is directly written to bar
```

```bash
echo nice > foo bar tar
```

```bash
cat foo
#nice bar tar
```

```bash
#except the file name, everything is written to the file foo
```

```bash
#Rule3, part2: whenever we encounter the |, the RHS is a command
```

```bash
date | tee foo bar
#Tue Jun 16 10:05:35 IST 2015
```

```bash
cat foo
#Tue Jun 16 10:05:35 IST 2015
```

```bash
cat bar
#Tue Jun 16 10:05:35 IST 2015
```

```bash
#command 'tee' give the output to a file/s and STDOUT too
```

```bash
cat /tmp/c.txt
#United Kingdom,57533000,London,6756000,Europe,1066,1945,-,English
#United States,252177000,Washington DC,606900,North America,1776,1945,-,English
#Venezuela,19733000,Caracas,1290087,Latin America,1811,1945,-,Spanish
#Vietnam,67589000,Hanoi,1088862,Orient,1945,1977,-,Vietnamese
#Yemen,1184300,San\'a,427150,Asia,1918,1957,Islam,Arabic
#Argentina,36955182,Buenos Aires,2033445,Latin America,1853,1945,-,Spanish
#Brazil,172860370,Brasilia,286037,Latin America,1822,1945,-,Portuguese
#Bahrain,634137,Manama,34137,Persian Gulf,1973,1977,Islamic,Arabic
#Cameroon,15421937,Yaounde,421937,Africa,1960,1974,-,Franch
#Djibouti,451442,Djibouti,1442,Africa,1977,1980,-,Franch
#Equatorial Guinea,474214,Malabo,74214,Africa,1991,1995,-,Franch
#Fiji,832494,Suva,32494,Oceania,1970,1975,-,English
#France,59329691,Paris,329691,Europe,486,1945,-,Franch
#Greece,10601527,Athens,601527,Europe,1829,1945,-,Greek
#Germany,82797408,Berlin,1797408,Europe,1871,1960,-,German
#Honduras,6249598,Tegucigalpa,1249598,Latin America,1821,1945,-,Spanish
#China,1261832482,Beijing,61832482,Asia,-221,1945,-,Chinese
#Canada,31281092,Ottawa,1281092,North America,1867,1945,-,English
#Hungary,10138844,Budapest,1138844,Europe,1001,1945,-,Hungerian
#India,1014003817,New Delhi,14003817,Asia,1947,1950,-,Indian
#Italy,57634327,Rome,3634327,Europe,1861,1950,-,Italian
#Ireland,3797257,Dublin,797257,Europe,1921,1945,-,English
#Japan,126549976,Tokio,16549976,Asia,-660,1955,-,Japanese
```

```bash
cat /etc/hosts
## Do not remove the following line, or various programs
## that require network functionality will fail.
#127.0.0.1               localhost.localdomain localhost qeserver.qetest
#::1             localhost6.localdomain6 localhost6 qeserver.qetest
#10.44.204.74 qeserver.qetest
```

```bash
cat -A /etc/hosts
## Do not remove the following line, or various programs$
## that require network functionality will fail.$
#127.0.0.1^I^Ilocalhost.localdomain localhost qeserver.qetest$
#::1^I^Ilocalhost6.localdomain6 localhost6 qeserver.qetest$
#10.44.204.74 qeserver.qetest$
```

```bash
#non-printable characters are shown with -A option
```

```bash
cat -n /etc/hosts
#     1  # Do not remove the following line, or various programs
#     2  # that require network functionality will fail.
#     3  127.0.0.1               localhost.localdomain localhost qeserver.qetest
#     4  ::1             localhost6.localdomain6 localhost6 qeserver.qetest
#     5  10.44.204.74 qeserver.qetest
```

```bash
#-n option shows the line numbers
```

```bash
head c.txt
#United Kingdom,57533000,London,6756000,Europe,1066,1945,-,English
#United States,252177000,Washington DC,606900,North America,1776,1945,-,English
#Venezuela,19733000,Caracas,1290087,Latin America,1811,1945,-,Spanish
#Vietnam,67589000,Hanoi,1088862,Orient,1945,1977,-,Vietnamese
#Yemen,1184300,San\'a,427150,Asia,1918,1957,Islam,Arabic
#Argentina,36955182,Buenos Aires,2033445,Latin America,1853,1945,-,Spanish
#Brazil,172860370,Brasilia,286037,Latin America,1822,1945,-,Portuguese
#Bahrain,634137,Manama,34137,Persian Gulf,1973,1977,Islamic,Arabic
#Cameroon,15421937,Yaounde,421937,Africa,1960,1974,-,Franch
#Djibouti,451442,Djibouti,1442,Africa,1977,1980,-,Franch
```

```bash
tail c.txt
#Greece,10601527,Athens,601527,Europe,1829,1945,-,Greek
#Germany,82797408,Berlin,1797408,Europe,1871,1960,-,German
#Honduras,6249598,Tegucigalpa,1249598,Latin America,1821,1945,-,Spanish
#China,1261832482,Beijing,61832482,Asia,-221,1945,-,Chinese
#Canada,31281092,Ottawa,1281092,North America,1867,1945,-,English
#Hungary,10138844,Budapest,1138844,Europe,1001,1945,-,Hungerian
#India,1014003817,New Delhi,14003817,Asia,1947,1950,-,Indian
#Italy,57634327,Rome,3634327,Europe,1861,1950,-,Italian
#Ireland,3797257,Dublin,797257,Europe,1921,1945,-,English
#Japan,126549976,Tokio,16549976,Asia,-660,1955,-,Japanese
```

```bash
head -5 c.txt
#United Kingdom,57533000,London,6756000,Europe,1066,1945,-,English
#United States,252177000,Washington DC,606900,North America,1776,1945,-,English
#Venezuela,19733000,Caracas,1290087,Latin America,1811,1945,-,Spanish
#Vietnam,67589000,Hanoi,1088862,Orient,1945,1977,-,Vietnamese
#Yemen,1184300,San\'a,427150,Asia,1918,1957,Islam,Arabic
```

```bash
tail -5 c.txt
#Hungary,10138844,Budapest,1138844,Europe,1001,1945,-,Hungerian
#India,1014003817,New Delhi,14003817,Asia,1947,1950,-,Indian
#Italy,57634327,Rome,3634327,Europe,1861,1950,-,Italian
#Ireland,3797257,Dublin,797257,Europe,1921,1945,-,English
#Japan,126549976,Tokio,16549976,Asia,-660,1955,-,Japanese
```


```bash
tail -n +20 c.txt
#India,1014003817,New Delhi,14003817,Asia,1947,1950,-,Indian
#Italy,57634327,Rome,3634327,Europe,1861,1950,-,Italian
#Ireland,3797257,Dublin,797257,Europe,1921,1945,-,English
#Japan,126549976,Tokio,16549976,Asia,-660,1955,-,Japanese
```

```bash
# '-n +20' will show all lines from line number 20
```

```bash
cut -d ',' -f1 c.txt #this shows the first column
#United Kingdom
#United States
#Venezuela
#Vietnam
#Yemen
#Argentina
#Brazil
#Bahrain
#Cameroon
#Djibouti
#Equatorial Guinea
#Fiji
#France
#Greece
#Germany
#Honduras
#China
#Canada
#Hungary
#India
#Italy
#Ireland
#Japan
```

```bash
cut -d ',' -f1,3 c.txt #this shows the 1st and 3rd column
#United Kingdom,London
#United States,Washington DC
#Venezuela,Caracas
#Vietnam,Hanoi
#Yemen,San\'a
#Argentina,Buenos Aires
#Brazil,Brasilia
#Bahrain,Manama
#Cameroon,Yaounde
#Djibouti,Djibouti
#Equatorial Guinea,Malabo
#Fiji,Suva
#France,Paris
#Greece,Athens
#Germany,Berlin
#Honduras,Tegucigalpa
#China,Beijing
#Canada,Ottawa
#Hungary,Budapest
#India,New Delhi
#Italy,Rome
#Ireland,Dublin
#Japan,Tokio
```

```bash
cut -d ',' -f1-3 c.txt #this shows the columns from 1-3
#United Kingdom,57533000,London
#United States,252177000,Washington DC
#Venezuela,19733000,Caracas
#Vietnam,67589000,Hanoi
#Yemen,1184300,San\'a
#Argentina,36955182,Buenos Aires
#Brazil,172860370,Brasilia
#Bahrain,634137,Manama
#Cameroon,15421937,Yaounde
#Djibouti,451442,Djibouti
#Equatorial Guinea,474214,Malabo
#Fiji,832494,Suva
#France,59329691,Paris
#Greece,10601527,Athens
#Germany,82797408,Berlin
#Honduras,6249598,Tegucigalpa
#China,1261832482,Beijing
#Canada,31281092,Ottawa
#Hungary,10138844,Budapest
#India,1014003817,New Delhi
#Italy,57634327,Rome
#Ireland,3797257,Dublin
#Japan,126549976,Tokio
```

```bash
cut -d ',' -f3- c.txt #this shows all the columns from 3
#London,6756000,Europe,1066,1945,-,English
#Washington DC,606900,North America,1776,1945,-,English
#Caracas,1290087,Latin America,1811,1945,-,Spanish
#Hanoi,1088862,Orient,1945,1977,-,Vietnamese
#San\'a,427150,Asia,1918,1957,Islam,Arabic
#Buenos Aires,2033445,Latin America,1853,1945,-,Spanish
#Brasilia,286037,Latin America,1822,1945,-,Portuguese
#Manama,34137,Persian Gulf,1973,1977,Islamic,Arabic
#Yaounde,421937,Africa,1960,1974,-,Franch
#Djibouti,1442,Africa,1977,1980,-,Franch
#Malabo,74214,Africa,1991,1995,-,Franch
#Suva,32494,Oceania,1970,1975,-,English
#Paris,329691,Europe,486,1945,-,Franch
#Athens,601527,Europe,1829,1945,-,Greek
#Berlin,1797408,Europe,1871,1960,-,German
#Tegucigalpa,1249598,Latin America,1821,1945,-,Spanish
#Beijing,61832482,Asia,-221,1945,-,Chinese
#Ottawa,1281092,North America,1867,1945,-,English
#Budapest,1138844,Europe,1001,1945,-,Hungerian
#New Delhi,14003817,Asia,1947,1950,-,Indian
#Rome,3634327,Europe,1861,1950,-,Italian
#Dublin,797257,Europe,1921,1945,-,English
#Tokio,16549976,Asia,-660,1955,-,Japanese
```

```bash
#non printable chars like TABS, NEW LINES are typed like this: Ctrl+v, TAB
```

```bash
echo 'a      b       c'
#a       b       c
```

```bash
echo 'a      b^Mc'
#c       b
```

### the tr command

```bash
cat c.txt | tr ',' '#'
#United Kingdom#57533000#London#6756000#Europe#1066#1945#-#English
#United States#252177000#Washington DC#606900#North America#1776#1945#-#English
#Venezuela#19733000#Caracas#1290087#Latin America#1811#1945#-#Spanish
#Vietnam#67589000#Hanoi#1088862#Orient#1945#1977#-#Vietnamese
#Yemen#1184300#San'a#427150#Asia#1918#1957#Islam#Arabic
#Argentina#36955182#Buenos Aires#2033445#Latin America#1853#1945#-#Spanish
#Brazil#172860370#Brasilia#286037#Latin America#1822#1945#-#Portuguese
#Bahrain#634137#Manama#34137#Persian Gulf#1973#1977#Islamic#Arabic
#Cameroon#15421937#Yaounde#421937#Africa#1960#1974#-#Franch
#Djibouti#451442#Djibouti#1442#Africa#1977#1980#-#Franch
#Equatorial Guinea#474214#Malabo#74214#Africa#1991#1995#-#Franch
#Fiji#832494#Suva#32494#Oceania#1970#1975#-#English
#France#59329691#Paris#329691#Europe#486#1945#-#Franch
#Greece#10601527#Athens#601527#Europe#1829#1945#-#Greek
#Germany#82797408#Berlin#1797408#Europe#1871#1960#-#German
#Honduras#6249598#Tegucigalpa#1249598#Latin America#1821#1945#-#Spanish
#China#1261832482#Beijing#61832482#Asia#-221#1945#-#Chinese
#Canada#31281092#Ottawa#1281092#North America#1867#1945#-#English
#Hungary#10138844#Budapest#1138844#Europe#1001#1945#-#Hungerian
#India#1014003817#New Delhi#14003817#Asia#1947#1950#-#Indian
#Italy#57634327#Rome#3634327#Europe#1861#1950#-#Italian
#Ireland#3797257#Dublin#797257#Europe#1921#1945#-#English
#Japan#126549976#Tokio#16549976#Asia#-660#1955#-#Japanese
```

```bash
cat c.txt | tr 'a-z' 'A-Z'
#UNITED KINGDOM,57533000,LONDON,6756000,EUROPE,1066,1945,-,ENGLISH
#UNITED STATES,252177000,WASHINGTON DC,606900,NORTH AMERICA,1776,1945,-,ENGLISH
#VENEZUELA,19733000,CARACAS,1290087,LATIN AMERICA,1811,1945,-,SPANISH
#VIETNAM,67589000,HANOI,1088862,ORIENT,1945,1977,-,VIETNAMESE
#YEMEN,1184300,SAN\'A,427150,ASIA,1918,1957,ISLAM,ARABIC
#ARGENTINA,36955182,BUENOS AIRES,2033445,LATIN AMERICA,1853,1945,-,SPANISH
#BRAZIL,172860370,BRASILIA,286037,LATIN AMERICA,1822,1945,-,PORTUGUESE
#BAHRAIN,634137,MANAMA,34137,PERSIAN GULF,1973,1977,ISLAMIC,ARABIC
#CAMEROON,15421937,YAOUNDE,421937,AFRICA,1960,1974,-,FRANCH
#DJIBOUTI,451442,DJIBOUTI,1442,AFRICA,1977,1980,-,FRANCH
#EQUATORIAL GUINEA,474214,MALABO,74214,AFRICA,1991,1995,-,FRANCH
#FIJI,832494,SUVA,32494,OCEANIA,1970,1975,-,ENGLISH
#FRANCE,59329691,PARIS,329691,EUROPE,486,1945,-,FRANCH
#GREECE,10601527,ATHENS,601527,EUROPE,1829,1945,-,GREEK
#GERMANY,82797408,BERLIN,1797408,EUROPE,1871,1960,-,GERMAN
#HONDURAS,6249598,TEGUCIGALPA,1249598,LATIN AMERICA,1821,1945,-,SPANISH
#CHINA,1261832482,BEIJING,61832482,ASIA,-221,1945,-,CHINESE
#CANADA,31281092,OTTAWA,1281092,NORTH AMERICA,1867,1945,-,ENGLISH
#HUNGARY,10138844,BUDAPEST,1138844,EUROPE,1001,1945,-,HUNGERIAN
#INDIA,1014003817,NEW DELHI,14003817,ASIA,1947,1950,-,INDIAN
#ITALY,57634327,ROME,3634327,EUROPE,1861,1950,-,ITALIAN
#IRELAND,3797257,DUBLIN,797257,EUROPE,1921,1945,-,ENGLISH
#JAPAN,126549976,TOKIO,16549976,ASIA,-660,1955,-,JAPANESE
```

```bash
echo hiiiiihiiiii | tr -s 'i'        #remove all occurances of 'i' but not one
#hihi
```

```bash
echo hiiiiihiiiii | tr -d 'i'        #remove all 'i'
#hh
```

```bash
cat -n c.txt | tail -n +10 | head -n 6
#    10  Djibouti,451442,Djibouti,1442,Africa,1977,1980,-,Franch
#    11  Equatorial Guinea,474214,Malabo,74214,Africa,1991,1995,-,Franch
#    12  Fiji,832494,Suva,32494,Oceania,1970,1975,-,English
#    13  France,59329691,Paris,329691,Europe,486,1945,-,Franch
#    14  Greece,10601527,Athens,601527,Europe,1829,1945,-,Greek
#    15  Germany,82797408,Berlin,1797408,Europe,1871,1960,-,German
```

```bash
cat -n c.txt | tail -n +10 | head -n 6
#    10  Djibouti,451442,Djibouti,1442,Africa,1977,1980,-,Franch
#    11  Equatorial Guinea,474214,Malabo,74214,Africa,1991,1995,-,Franch
#    12  Fiji,832494,Suva,32494,Oceania,1970,1975,-,English
#    13  France,59329691,Paris,329691,Europe,486,1945,-,Franch
#    14  Greece,10601527,Athens,601527,Europe,1829,1945,-,Greek
#    15  Germany,82797408,Berlin,1797408,Europe,1871,1960,-,German
```

```bash
#this is not more efficient as we are reading the whole file.
```

```bash
head -15 c.txt | cat -n | tail -6    #more efficient command
#    10  Djibouti,451442,Djibouti,1442,Africa,1977,1980,-,Franch
#    11  Equatorial Guinea,474214,Malabo,74214,Africa,1991,1995,-,Franch
#    12  Fiji,832494,Suva,32494,Oceania,1970,1975,-,English
#    13  France,59329691,Paris,329691,Europe,486,1945,-,Franch
#    14  Greece,10601527,Athens,601527,Europe,1829,1945,-,Greek
#    15  Germany,82797408,Berlin,1797408,Europe,1871,1960,-,German
```

```bash
cal
#     June 2015
#Su Mo Tu We Th Fr Sa
#    1  2  3  4  5  6
# 7  8  9 10 11 12 13
#14 15 16 17 18 19 20
#21 22 23 24 25 26 27
#28 29 30
#
```

```bash
cal | head -5
#     June 2015
#Su Mo Tu We Th Fr Sa
#    1  2  3  4  5  6
# 7  8  9 10 11 12 13
#14 15 16 17 18 19 20

cal | head -5 | tail -1 | cut -d ' ' -f2
#15

#show 15 from the cal

cal | head -5 | tail -2 | cut -d ' ' -f2
#7
#15

cal | head -6 | tail -2 | cut -d ' ' -f2
#15
#22

cal | head -6 | tail -2 | cut -d ' ' -f2 | tr '\n' ' '
#15 22

#now there is no new line after 22
#so we use command substutition

echo `cal | head -6 | tail -2 | cut -d ' ' -f2 | tr '\n' ' '`
#15 22

#we cannot use echo as it takes input

cal | head -6 | tail -2 | cut -d ' ' -f2 | tr '\n' ' ' | echo
#

#piping it to echo did not gave the output
#we use xargs which gives input to echo

cal | head -6 | tail -2 | cut -d ' ' -f2 | tr '\n' ' ' | xargs echo
#15 22
```

```bash
df
#Filesystem           1K-blocks      Used Available Use% Mounted on
#/dev/sda1            147386176  16589920 123188680  12% /
#tmpfs                  1032324         0   1032324   0% /dev/shm


#task: show 12 in this command output

df | tr -s ' '
#Filesystem 1K-blocks Used Available Use% Mounted on
#/dev/sda1 147386176 16589936 123188664 12% /
#tmpfs 1032324 0 1032324 0% /dev/shm

df | head -2 | tr -s ' ' |
#>

df | head -2 | tr -s ' '
#Filesystem 1K-blocks Used Available Use% Mounted on
#/dev/sda1 147386176 16589944 123188656 12% /

df | head -2 | tail -1 | tr -s ' '
#/dev/sda1 147386176 16589944 123188656 12% /

df | head -2 | tail -1 | tr -s ' ' | cut -d ' ' -f5
#12%

df | head -2 | tail -1 | tr -s ' ' | cut -d ' ' -f5 | cut -d '%' -f1
#12

df | head -2 | tail -1 | tr -s ' ' | cut -d ' ' -f5 | tr -d '%'
#12
```

```bash
/sbin/ifconfig eth1
#eth1      Link encap:Ethernet  HWaddr 00:1B:11:15:E5:D8
#          inet addr:10.44.204.213  Bcast:10.44.207.255  Mask:255.255.252.0
#          inet6 addr: fe80::21b:11ff:fe15:e5d8/64 Scope:Link
#          UP BROADCAST RUNNING MULTICAST  MTU:1500  Metric:1
#          RX packets:4601250 errors:0 dropped:0 overruns:0 frame:0
#          TX packets:359836 errors:0 dropped:0 overruns:0 carrier:0
#          collisions:0 txqueuelen:1000
#          RX bytes:357173011 (340.6 MiB)  TX bytes:42759705 (40.7 MiB)
#          Interrupt:217 Base address:0xef00
#


#task: show the ip address after 'inet addr:'

/sbin/ifconfig eth1 | head -2 | tail -1
#          inet addr:10.44.204.213  Bcast:10.44.207.255  Mask:255.255.252.0

/sbin/ifconfig eth1 | head -2 | tail -1 | tr -s ' '
# inet addr:10.44.204.213 Bcast:10.44.207.255 Mask:255.255.252.0

/sbin/ifconfig eth1 | head -2 | tail -1 | tr -s ' ' | cut -d ' ' -f2
#inet

/sbin/ifconfig eth1 | head -2 | tail -1 | tr -s ' ' | cut -d ' ' -f3
#addr:10.44.204.213

/sbin/ifconfig eth1 | head -2 | tail -1 | tr -s ' ' | cut -d ' ' -f3 | cut -d ':' -f2
#10.44.204.213

/sbin/ifconfig eth1 | head -2 | tail -1 | cut -d ':' -f2 | cut -d ' ' -f1
#10.44.204.213
```

## Grep command

```bash
grep 'India' c.txt
#India,1014003817,New Delhi,14003817,Asia,1947,1950,-,Indian
```

```bash
grep 'English' c.txt
#United Kingdom,57533000,London,6756000,Europe,1066,1945,-,English
#United States,252177000,Washington DC,606900,North America,1776,1945,-,English
#Fiji,832494,Suva,32494,Oceania,1970,1975,-,English
#Canada,31281092,Ottawa,1281092,North America,1867,1945,-,English
#Ireland,3797257,Dublin,797257,Europe,1921,1945,-,English
```

```bash
grep --color -i 'english' c.txt
#United Kingdom,57533000,London,6756000,Europe,1066,1945,-,English
#United States,252177000,Washington DC,606900,North America,1776,1945,-,English
#Fiji,832494,Suva,32494,Oceania,1970,1975,-,English
#Canada,31281092,Ottawa,1281092,North America,1867,1945,-,English
#Ireland,3797257,Dublin,797257,Europe,1921,1945,-,English
```

```bash
#show all lines which do not contain 'english'
grep --color -i -v 'english' c.txt   
#Venezuela,19733000,Caracas,1290087,Latin America,1811,1945,-,Spanish
#Vietnam,67589000,Hanoi,1088862,Orient,1945,1977,-,Vietnamese
#Yemen,1184300,San\'a,427150,Asia,1918,1957,Islam,Arabic
#Argentina,36955182,Buenos Aires,2033445,Latin America,1853,1945,-,Spanish
#Brazil,172860370,Brasilia,286037,Latin America,1822,1945,-,Portuguese
#Bahrain,634137,Manama,34137,Persian Gulf,1973,1977,Islamic,Arabic
#Cameroon,15421937,Yaounde,421937,Africa,1960,1974,-,Franch
#Djibouti,451442,Djibouti,1442,Africa,1977,1980,-,Franch
#Equatorial Guinea,474214,Malabo,74214,Africa,1991,1995,-,Franch
#France,59329691,Paris,329691,Europe,486,1945,-,Franch
#Greece,10601527,Athens,601527,Europe,1829,1945,-,Greek
#Germany,82797408,Berlin,1797408,Europe,1871,1960,-,German
#Honduras,6249598,Tegucigalpa,1249598,Latin America,1821,1945,-,Spanish
#China,1261832482,Beijing,61832482,Asia,-221,1945,-,Chinese
#Hungary,10138844,Budapest,1138844,Europe,1001,1945,-,Hungerian
#India,1014003817,New Delhi,14003817,Asia,1947,1950,-,Indian
#Italy,57634327,Rome,3634327,Europe,1861,1950,-,Italian
#Japan,126549976,Tokio,16549976,Asia,-660,1955,-,Japanese
```

```bash
#count all lines which do not contain 'english'
grep --color -i -v -c 'english' c.txt        
#18
```

```bash
#show line containing India and one line after that
grep --color -A1 'India' c.txt       
#India,1014003817,New Delhi,14003817,Asia,1947,1950,-,Indian
#Italy,57634327,Rome,3634327,Europe,1861,1950,-,Italian
```

```bash
#show line containing India and one line before that
grep --color -B1 'India' c.txt       
#Hungary,10138844,Budapest,1138844,Europe,1001,1945,-,Hungerian
#India,1014003817,New Delhi,14003817,Asia,1947,1950,-,Indian
```

```bash
#show line containing India and one line before and after that
grep --color -C1 'India' c.txt       
#Hungary,10138844,Budapest,1138844,Europe,1001,1945,-,Hungerian
#India,1014003817,New Delhi,14003817,Asia,1947,1950,-,Indian
#Italy,57634327,Rome,3634327,Europe,1861,1950,-,Italian
```

```bash
egrep 'English|Spanish' c.txt
#United Kingdom,57533000,London,6756000,Europe,1066,1945,-,English
#United States,252177000,Washington DC,606900,North America,1776,1945,-,English
#Venezuela,19733000,Caracas,1290087,Latin America,1811,1945,-,Spanish
#Argentina,36955182,Buenos Aires,2033445,Latin America,1853,1945,-,Spanish
#Fiji,832494,Suva,32494,Oceania,1970,1975,-,English
#Honduras,6249598,Tegucigalpa,1249598,Latin America,1821,1945,-,Spanish
#Canada,31281092,Ottawa,1281092,North America,1867,1945,-,English
#Ireland,3797257,Dublin,797257,Europe,1921,1945,-,English
```

```bash
egrep '(Engl|Span)ish' c.txt
#United Kingdom,57533000,London,6756000,Europe,1066,1945,-,English
#United States,252177000,Washington DC,606900,North America,1776,1945,-,English
#Venezuela,19733000,Caracas,1290087,Latin America,1811,1945,-,Spanish
#Argentina,36955182,Buenos Aires,2033445,Latin America,1853,1945,-,Spanish
#Fiji,832494,Suva,32494,Oceania,1970,1975,-,English
#Honduras,6249598,Tegucigalpa,1249598,Latin America,1821,1945,-,Spanish
#Canada,31281092,Ottawa,1281092,North America,1867,1945,-,English
#Ireland,3797257,Dublin,797257,Europe,1921,1945,-,English
```

```bash
egrep '^I' c.txt
#India,1014003817,New Delhi,14003817,Asia,1947,1950,-,Indian
#Italy,57634327,Rome,3634327,Europe,1861,1950,-,Italian
#Ireland,3797257,Dublin,797257,Europe,1921,1945,-,English
```

```bash
egrep 'ian$' c.txt
#Hungary,10138844,Budapest,1138844,Europe,1001,1945,-,Hungerian
#India,1014003817,New Delhi,14003817,Asia,1947,1950,-,Indian
#Italy,57634327,Rome,3634327,Europe,1861,1950,-,Italian
```

```bash
egrep '^I[nt]' c.txt
#India,1014003817,New Delhi,14003817,Asia,1947,1950,-,Indian
#Italy,57634327,Rome,3634327,Europe,1861,1950,-,Italian
```

```bash
egrep '^I[^nt]' c.txt
#Ireland,3797257,Dublin,797257,Europe,1921,1945,-,English
```

```bash
egrep '^I[n-t]' c.txt
#India,1014003817,New Delhi,14003817,Asia,1947,1950,-,Indian
#Italy,57634327,Rome,3634327,Europe,1861,1950,-,Italian
#Ireland,3797257,Dublin,797257,Europe,1921,1945,-,English
```

```bash
egrep '^I....' c.txt
#India,1014003817,New Delhi,14003817,Asia,1947,1950,-,Indian
#Italy,57634327,Rome,3634327,Europe,1861,1950,-,Italian
#Ireland,3797257,Dublin,797257,Europe,1921,1945,-,English
```

```bash
egrep '^I.{4}' c.txt
#India,1014003817,New Delhi,14003817,Asia,1947,1950,-,Indian
#Italy,57634327,Rome,3634327,Europe,1861,1950,-,Italian
#Ireland,3797257,Dublin,797257,Europe,1921,1945,-,English
```

```bash
egrep '^I.{4},' c.txt
#India,1014003817,New Delhi,14003817,Asia,1947,1950,-,Indian
#Italy,57634327,Rome,3634327,Europe,1861,1950,-,Italian
```

```bash
egrep '^I.{4,6},' c.txt
#India,1014003817,New Delhi,14003817,Asia,1947,1950,-,Indian
#Italy,57634327,Rome,3634327,Europe,1861,1950,-,Italian
#Ireland,3797257,Dublin,797257,Europe,1921,1945,-,English
```

```bash
alias egrep='egrep --color'
```

```bash
egrep '^.{5},[0-9]{10},' c.txt
#China,1261832482,Beijing,61832482,Asia,-221,1945,-,Chinese
#India,1014003817,New Delhi,14003817,Asia,1947,1950,-,Indian
```

```bash
egrep '^ab{0,1}c$'
#ac
#ac
#abc
#abc
#abbc
```

```bash
egrep '^ab?c$'
#ac
#ac
#abc
#abc
#ab
```

```bash
egrep '^ab+c$'
#accb
#abbbc
#abbbc
#ac
#ac
```

```bash
egrep 'gr[ae]y colou?r'
#gray color
#gray color
#grey color
#grey color
#gey cik
```

```bash
egrep '(wo)?m[ae]n of honou?r'
#men of honor
#men of honor
#man og honour
#man of honor
#man of honor
```

```bash
egrep '^[gj]eo?f{2}(er|re)y'
#geoffery
#geoffery
#geoffrey
#geoffrey
#jeffery
#jeffery
#jeffrey
#jeffrey
#
#
#jeoffery
#jeoffery
```

```bash
egrep '^(geo|je)f{2}(er|re)y$'
#jeoffery
#geoffery
#geoffery
```

```bash
cat fox-tree.dat
#The brown fox climbed the brown tree.
#The brown fox climbed the grey tree.
#The brown fox climbed the brown tree.
#The grey fox climbed the brown tree.
#The grey fox climbed the grey tree.
#The yellow fox climbed the yellow tree.
#The yellow fox climbed the red tree.
#The red fox climbed the red tree.
#The red fox climbed the brown tree.
#The red fox climbed the grey tree.
```

```bash
egrep 'The (.*) fox.*\1 tree' fox-tree.dat
#The brown fox climbed the brown tree.
#The brown fox climbed the brown tree.
#The grey fox climbed the grey tree.
#The yellow fox climbed the yellow tree.
#The red fox climbed the red tree.
```


## find PATHs CRITERIAs


```bash
find ./ -name '*.txt'
#./c.txt
#./file.txt
#./tmp.txt
```

```bash
find -type d
#.
#./.mozilla
#./.mozilla/extensions
#./.mozilla/plugins
#./.kde
#./.kde/Autostart
```

```bash
find ./ -type d
#./
#./.mozilla
#./.mozilla/extensions
#./.mozilla/plugins
#./.kde
#./.kde/Autostart
```

```bash
find ./ -type f
#./.bashrc
#./bar
#./c.txt
#./.emacs
#./file.txt
#./tmp.txt
#./.bash_logout
#./fox-tree.dat
#./.bash_history
#./foo
#./.bash_profile
#./.zshrc
#./.kde/Autostart/.directory
#./H
#./.viminfo
#./fun_file
```

```bash
find ./ ! -type f
#./
#./.mozilla
#./.mozilla/extensions
#./.mozilla/plugins
#./.kde
#./.kde/Autostart
```

```bash
find ./ -size +100c
#./
#./.bashrc
#./c.txt
#./.mozilla
#./.mozilla/extensions
#./.mozilla/plugins
#./.emacs
#./file.txt
#./fox-tree.dat
#./.bash_history
#./.bash_profile
#./.zshrc
#./.kde
#./.kde/Autostart
#./.kde/Autostart/.directory
#./.viminfo
```

## Shell scripting

```bash
vim first.bash

cat first.bash
#echo 'hello world'
```

```bash
bash first.bash
#hello world
```

```bash
mv first.bash first

cat first
#echo 'hello world'
```

```bash
vim first

cat first
##!/bin/ksh
#
#echo 'hello world'
```

```bash
chmod +x first

./first
#hello world
```

```bash
vim hello-read.bash

cat hello-read.bash
##!/bin/bash
#
#echo 'Please enter a name'
#read name
#echo "Hello ${name}"
```

```bash
ls
#bar  c.txt  file.txt  first  foo  fox-tree.dat  fun_file  H  hello-read.bash  sample1  tmp.txt
```

```bash
bash ./hello-read.bash
#Please enter a name
#lakhan
#Hello lakhan
```

```bash
vim hello-args.bash

cat hello-args.bash
##!/bin/bash
#
#echo "hello $1"
#echo "bye $2"
```

```bash
./hello-args.bash aaa bbb
#hello aaa
#bye bbb
```

```bash
vim args.bash

cat args.bash
##!/bin/bash
#
#echo "$1---$2---$3"
#echo "Total -> $#"
#echo "All -> $*"
#echo "Name of the script -> $0"
```

```bash
chmod +x ./args.bash

./args.bash aaa bbb ccc
#aaa---bbb---ccc
#Total -> 3
#All -> aaa bbb ccc
#Name of the script -> ./args.bash
```

#### rule4: always give a command after 'if'

```bash
vim if0.bash

cat if0.bash
##!/bin/bash
#
#if date
#then
#        echo yes
#else
#        echo no
#fi
```

```bash
bash ./if0.bash
#Tue Jun 16 14:55:46 IST 2015
#yes
```

```bash
vim checkUser.bash

cat checkUser.bash
##!/bin/bash
#
#echo "Who are you?"
#read name
#
#if who | egrep "^${name}" > /dev/null
#then
#        echo "You are ${name} and logged in"
#else
#        echo "You are ${name} and not logged in currently"
#fi
```

```bash
bash checkUser.bash
#Who are you?
#bashshell14
#You are bashshell14 and logged in
```

```bash
vim numericCompare.bash

cat numericCompare.bash
##!/bin/bash
#
#if test $# -ne 2
#then
#        echo 'Incorrect usage'
#        exit 1
#fi
#
#
#if test ${1} -gt ${2}
#then
#        echo "${1} > ${2}"
#elif test ${1} -eq ${2}
#then
#        echo "${1} = ${2}"
#else
#        echo "${1} < ${2}"
#fi
```

```bash
./numericCompare.bash 4 5
#4 < 5

./numericCompare.bash 4 4
#4 = 4

./numericCompare.bash 4 1
#4 > 1
```

```bash
vim numericCompare2.bash

cat numericCompare2.bash
##!/bin/bash
#
#if [ $# -ne 2 ]
#then
#        echo 'Incorrect usage'
#        exit 1
#fi
#
#
#if [ ${1} -gt ${2} ]
#then
#        echo "${1} > ${2}"
#elif [ ${1} -eq ${2} ]
#then
#        echo "${1} = ${2}"
#else
#        echo "${1} < ${2}"
#fi
```

```bash
chmod +x numericCompare2.bash

./numericCompare2.bash 5 6
#5 < 6
```

```bash
vim stringComp.bash
chmod +x stringComp.bash

cat stringComp.bash
##!/bin/bash
#
#echo 'enter the name'
#read name
#
#if [[ -z $name ]]
#then
#        echo 'nothing entered'
#        exit 1
#fi
#
#name=`echo "$name" | tr 'A-Z' 'a-z'`
#
#if [[ $name == lakhan ]]
#then
#        echo "$name == lakhan"
#elif [[ "$name" > lakhan ]]
#then
#        echo  "$name > lakhan"
#else
#        echo  "$name < lakhan"
#fi
```

```bash
./stringComp.bash
#enter the name
#lakhan
#lakhan == lakhan

./stringComp.bash
#enter the name
#test
#test > lakhan

./stringComp.bash
#enter the name
#abc
#abc < lakhan
```

#### Rule5.1: [	] and [[	]] are same except: 1) quoting is not compulsory in the [[	]]


```bash
mkfifo /tmp/lakhan.fifo

ls -l /tmp/lakhan.fifo
#prw-rw-r-- 1 bashshell14 bashshell14 0 Jun 16 16:25 /tmp/lakhan.fifo

chmod 777 /tmp/lakhan.fifo

cat /etc/hosts > /tmp/lakhan.fifo

#the other user will type following command

tr 'a-z' 'A-Z' < /tmp/lakhan.fifo
```

```bash
vim compareThreeNos.bash

chmod +x compareThreeNos.bash

cat compareThreeNos.bash
##!/bin/bash
#
#if [[ $1 -gt $2 ]] && [[ $2 -gt $3 ]]
##if [[  $1 -gt $2 && $2 -gt $3 ]]       : this is also allowed
##if [ $1 -gt $2 ] && [ $2 -gt $3 ]      : this is also allowed
##if [ $1 -gt $2 && $2 -gt $3 ]          : this is not allowed
#then
#        echo yes
#else
#        echo no
#fi
```

```bash
./compareThreeNos.bash 12 13 14
#no

./compareThreeNos.bash 12 11 10
#yes
```

## sed 'FILTER' filename

```bash
#sed by default takes input and gives output, similar to cat
#it conditionally applies filter. if filter is not applied, it works similar to cat
#filter has two parts
#sed is a line operator
#it reads file line by line
#it takes the first line of file and stores in the pattern space
#then it performes substitution
#then it prints it to the screen
#erase the pattern space to take next line
#if the pattern is not found in line, it prints it as it is
```

```bash
cat d.txt
#In Xanadu did Kubla Khan
#A stately pleasure dome decree:
#Where Alph, the sacred river, ran
#Through caverns measureless to man
#Down to a sunless sea.
```

```bash
sed 's/ /#/g' d.txt
#In#Xanadu#did#Kubla#Khan
#A#stately#pleasure#dome#decree:
#Where#Alph,#the#sacred#river,#ran
#Through#caverns#measureless#to#man
#Down#to#a#sunless#sea.
```

```bash
#replaces 2nd space with #
sed 's/ /#/2' d.txt          
#In Xanadu#did Kubla Khan
#A stately#pleasure dome decree:
#Where Alph,#the sacred river, ran
#Through caverns#measureless to man
#Down to#a sunless sea.
```

```bash
#replaces all spaces after 3rd with #
sed 's/ /#/3g' d.txt                 
#In Xanadu did#Kubla#Khan
#A stately pleasure#dome#decree:
#Where Alph, the#sacred#river,#ran
#Through caverns measureless#to#man
#Down to a#sunless#sea.
```

```bash
#ignore case
sed 's/khan/Don/i' d.txt             
#In Xanadu did Kubla Don
#A stately pleasure dome decree:
#Where Alph, the sacred river, ran
#Through caverns measureless to man
#Down to a sunless sea.
```

```bash
cat e.txt
#SR      8649    275     Asia
#Canada  3852    25      North America
#China   3705    1032    Asia
#USA     3615    237     North America
#Brazil  3286    134     South America
#India   1267    746     Asia
#Mexico  762     78      North America
#France  211     55      Europe
#Japan   144     120     Asia
#Germany 96      61      Europe
#England 94      56      Europe
```

```bash
sed 's/[0-9]/-/g' e.txt
#SR      ----    ---     Asia
#Canada  ----    --      North America
#China   ----    ----    Asia
#USA     ----    ---     North America
#Brazil  ----    ---     South America
#India   ----    ---     Asia
#Mexico  ---     --      North America
#France  ---     --      Europe
#Japan   ---     ---     Asia
#Germany --      --      Europe
#England --      --      Europe
```

```bash
#perform this action on 2nd line only
sed '2s/[0-9]/-/g' e.txt     
#SR      8649    275     Asia
#Canada  ----    --      North America
#China   3705    1032    Asia
#USA     3615    237     North America
#Brazil  3286    134     South America
#India   1267    746     Asia
#Mexico  762     78      North America
#France  211     55      Europe
#Japan   144     120     Asia
#Germany 96      61      Europe
#England 94      56      Europe
```

```bash
#search for a pattern and perform the action on that line only
sed '/Canada/ s/[0-9]/-/g' e.txt     
#SR      8649    275     Asia
#Canada  ----    --      North America
#China   3705    1032    Asia
#USA     3615    237     North America
#Brazil  3286    134     South America
#India   1267    746     Asia
#Mexico  762     78      North America
#France  211     55      Europe
#Japan   144     120     Asia
#Germany 96      61      Europe
#England 94      56      Europe
```

```bash
#search for a pattern and perform the action on that line only
sed '/Canada/ s/[0-9]/-/g' e.txt     
#SR      8649    275     Asia
#Canada  ----    --      North America
#China   3705    1032    Asia
#USA     3615    237     North America
#Brazil  3286    134     South America
#India   1267    746     Asia
#Mexico  762     78      North America
#France  211     55      Europe
#Japan   144     120     Asia
#Germany 96      61      Europe
#England 94      56      Europe
```

```bash
#perform this action from 2nd row to 5th row
sed '2,5 s/[0-9]/-/g' e.txt  
#SR      8649    275     Asia
#Canada  ----    --      North America
#China   ----    ----    Asia
#USA     ----    ---     North America
#Brazil  ----    ---     South America
#India   1267    746     Asia
#Mexico  762     78      North America
#France  211     55      Europe
#Japan   144     120     Asia
#Germany 96      61      Europe
#England 94      56      Europe
```

```bash
#search for a pattern and perform the action till the line number provided
sed '/Canada/,5 s/[0-9]/-/g' e.txt   
#SR      8649    275     Asia
#Canada  ----    --      North America
#China   ----    ----    Asia
#USA     ----    ---     North America
#Brazil  ----    ---     South America
#India   1267    746     Asia
#Mexico  762     78      North America
#France  211     55      Europe
#Japan   144     120     Asia
#Germany 96      61      Europe
#England 94      56      Europe
```

```bash
#search for a pattern and perform the action till next pattern found
sed '/North/,/South/ s/[0-9]/-/g' e.txt      
#SR      8649    275     Asia
#Canada  ----    --      North America
#China   ----    ----    Asia
#USA     ----    ---     North America
#Brazil  ----    ---     South America
#India   1267    746     Asia
#Mexico  ---     --      North America
#France  ---     --      Europe
#Japan   ---     ---     Asia
#Germany --      --      Europe
#England --      --      Europe
```

```bash
#perform this action from 2nd row till last row
sed '2,$ s/[0-9]/-/g' e.txt  
#SR      8649    275     Asia
#Canada  ----    --      North America
#China   ----    ----    Asia
#USA     ----    ---     North America
#Brazil  ----    ---     South America
#India   ----    ---     Asia
#Mexico  ---     --      North America
#France  ---     --      Europe
#Japan   ---     ---     Asia
#Germany --      --      Europe
#England --      --      Europe
```

```bash
sed '2! s/[0-9]/-/g' e.txt   #negation
#SR      ----    ---     Asia
#Canada  3852    25      North America
#China   ----    ----    Asia
#USA     ----    ---     North America
#Brazil  ----    ---     South America
#India   ----    ---     Asia
#Mexico  ---     --      North America
#France  ---     --      Europe
#Japan   ---     ---     Asia
#Germany --      --      Europe
#England --      --      Europe
```

```bash
#negation
sed '/North/,/South/! s/[0-9]/-/g' e.txt     
#SR      ----    ---     Asia
#Canada  3852    25      North America
#China   3705    1032    Asia
#USA     3615    237     North America
#Brazil  3286    134     South America
#India   ----    ---     Asia
#Mexico  762     78      North America
#France  211     55      Europe
#Japan   144     120     Asia
#Germany 96      61      Europe
#England 94      56      Europe
```

```bash
#this prints all the lines
sed 's/Khan/Don/' d.txt      
#In Xanadu did Kubla Don
#A stately pleasure dome decree:
#Where Alph, the sacred river, ran
#Through caverns measureless to man
#Down to a sunless sea.
```

```bash
#no lines are printed
sed -n 's/Khan/Don/' d.txt   
```

```bash
#print only the lines for which substitution was done
sed -n 's/Khan/Don/p' d.txt 
#In Xanadu did Kubla Don
```

```bash
grep 'Khan' d.txt
#In Xanadu did Kubla Khan
```

```bash
#use of sed similar to grep
sed -n '/Khan/p' d.txt       
#In Xanadu did Kubla Khan
```

```bash
grep -v 'Khan' d.txt
#A stately pleasure dome decree:
#Where Alph, the sacred river, ran
#Through caverns measureless to man
#Down to a sunless sea.
```

```bash
sed '/Khan/d' d.txt
#A stately pleasure dome decree:
#Where Alph, the sacred river, ran
#Through caverns measureless to man
#Down to a sunless sea.
```

```bash
sed -n '/Asia/p' e.txt
#SR      8649    275     Asia
#China   3705    1032    Asia
#India   1267    746     Asia
#Japan   144     120     Asia
```

```bash
#search for the line with pattern Canada and print line next to it
sed -n '/Canada/{n;p;}' e.txt        
#China   3705    1032    Asia
```

```bash
#search for the line with pattern Canada and print that line and next line
sed -n '/Canada/{p;n;p;}' e.txt      
#Canada  3852    25      North America
#China   3705    1032    Asia
```

```bash
#search for the line with pattern Canada and print that line and next line: N appends the second line to pattern space
sed -n '/Canada/{N;p;}' e.txt        
#Canada  3852    25      North America
#China   3705    1032    Asia
```

```bash
sed -n '/Asia/p' e.txt
#SR      8649    275     Asia
#China   3705    1032    Asia
#India   1267    746     Asia
#Japan   144     120     Asia
```

```bash
#search for the line with pattern Canada and print line next to it
sed -n '/Canada/{n;p;}' e.txt        
#China   3705    1032    Asia
```

```bash
#search for the line with pattern Canada and print that line and next line
sed -n '/Canada/{p;n;p;}' e.txt      
#Canada  3852    25      North America
#China   3705    1032    Asia
```

```bash
#search for the line with pattern Canada and print that line and next line: N appends the second line to pattern space
sed -n '/Canada/{N;p;}' e.txt        
#Canada  3852    25      North America
#China   3705    1032    Asia
```

```bash
#prints line from 2-5
sed -n '2,5p' e.txt          
#Canada  3852    25      North America
#China   3705    1032    Asia
#USA     3615    237     North America
#Brazil  3286    134     South America
```

```bash
#prints line #2 and line #5
sed -n -e '2p' -e '5p' e.txt                 
#Canada  3852    25      North America
#Brazil  3286    134     South America
```

```bash
#different type of substitutions from line 2 to line 5
sed -n -e '2s/[0-9]/-/gp' -e 's/ /#/gp' e.txt                
#Canada  ----    --      North America
#Canada  ----    --      North#America
#USA     3615    237     North#America
#Brazil  3286    134     South#America
#Mexico  762     78      North#America
```

```bash
#different type of substitutions on line 2 and line 5
sed -n -e '2s/[0-9]/-/gp' -e '5s/ /#/gp' e.txt               
#Canada  ----    --      North America
#Brazil  3286    134     South#America
```

```bash
echo a:b:c:d:e:f | sed -e 's/:/#/2' -e 's/:/#/3'
#a:b#c:d#e:f

#replace 2nd and 4th : with #
```

#### sed: backreferences

```bash
echo shantanu:kulkarni
#shantanu:kulkarni

echo shantanu:kulkarni | sed -r 's/shantanu:kulkarni/kulkarni:shantanu/'
#kulkarni:shantanu

echo shantanu:kulkarni | sed -r 's/(shantanu):(kulkarni)/kulkarni:shantanu/'
#kulkarni:shantanu

echo shantanu:kulkarni | sed -r 's/(shantanu):(kulkarni)/\1:\2/'
#shantanu:kulkarni

echo shantanu:kulkarni | sed -r 's/(........):(........)/\2:\1/'
#shantanu:kulkarni

echo shantram:kulkarni | sed -r 's/(........):(........)/\2:\1/'
#shantram:kulkarni

echo shantram:kulkarni | sed -r 's/(........):(........)/\2:\1/'
#kulkarni:shantram

echo shantram:kulkarni | sed -r 's/(.*):(.*)/\2:\1/'
#kulkarni:shantram

echo anshuman | sed -r 's/(ans)(human)/\1\2 is a \2/'
#anshuman is a human

echo 'pushkar nagarkar' | sed -r 's/(.*) (nagar)kar/\1 is not from \2/'
#pushkar is not from nagar

echo 'shantanu' | sed -r 's/(s.*a)nu/\1ram/'
#shantaram

echo 'shantanu' | sed -r 's/(s.*a)n/\1ram/'
#shantaramu

echo 'hello world' | sed 's/h.*o/hi/'
#hirld

echo hello hi ok globe | sed -r 's/(.*) (.*)/\1#\2/'
#hello hi ok#globe
```

```bash
cat c.txt
#United Kingdom,57533000,London,6756000,Europe,1066,1945,-,English
#United States,252177000,Washington DC,606900,North America,1776,1945,-,English
#Venezuela,19733000,Caracas,1290087,Latin America,1811,1945,-,Spanish
#Vietnam,67589000,Hanoi,1088862,Orient,1945,1977,-,Vietnamese
#Yemen,1184300,San\'a,427150,Asia,1918,1957,Islam,Arabic
#Argentina,36955182,Buenos Aires,2033445,Latin America,1853,1945,-,Spanish
#Brazil,172860370,Brasilia,286037,Latin America,1822,1945,-,Portuguese
#Bahrain,634137,Manama,34137,Persian Gulf,1973,1977,Islamic,Arabic
#Cameroon,15421937,Yaounde,421937,Africa,1960,1974,-,Franch
#Djibouti,451442,Djibouti,1442,Africa,1977,1980,-,Franch
#Equatorial Guinea,474214,Malabo,74214,Africa,1991,1995,-,Franch
#Fiji,832494,Suva,32494,Oceania,1970,1975,-,English
#France,59329691,Paris,329691,Europe,486,1945,-,Franch
#Greece,10601527,Athens,601527,Europe,1829,1945,-,Greek
#Germany,82797408,Berlin,1797408,Europe,1871,1960,-,German
#Honduras,6249598,Tegucigalpa,1249598,Latin America,1821,1945,-,Spanish
#China,1261832482,Beijing,61832482,Asia,-221,1945,-,Chinese
#Canada,31281092,Ottawa,1281092,North America,1867,1945,-,English
#Hungary,10138844,Budapest,1138844,Europe,1001,1945,-,Hungerian
#India,1014003817,New Delhi,14003817,Asia,1947,1950,-,Indian
#Italy,57634327,Rome,3634327,Europe,1861,1950,-,Italian
#Ireland,3797257,Dublin,797257,Europe,1921,1945,-,English
#Japan,126549976,Tokio,16549976,Asia,-660,1955,-,Japanese
```

```bash
sed -n '10,15p' c.txt
#Djibouti,451442,Djibouti,1442,Africa,1977,1980,-,Franch
#Equatorial Guinea,474214,Malabo,74214,Africa,1991,1995,-,Franch
#Fiji,832494,Suva,32494,Oceania,1970,1975,-,English
#France,59329691,Paris,329691,Europe,486,1945,-,Franch
#Greece,10601527,Athens,601527,Europe,1829,1945,-,Greek
#Germany,82797408,Berlin,1797408,Europe,1871,1960,-,German
```

```bash
cal
#     June 2015
#Su Mo Tu We Th Fr Sa
#    1  2  3  4  5  6
# 7  8  9 10 11 12 13
#14 15 16 17 18 19 20
#21 22 23 24 25 26 27
#28 29 30
#
```

```bash
cal | sed -n -r '5s/...(..).*/\1/p'
#15
```

```bash
#print 15 from the cal
```

```bash
df
#Filesystem           1K-blocks      Used Available Use% Mounted on
#/dev/sda1            147386176  16617400 123161200  12% /
#tmpfs                  1032324         0   1032324   0% /dev/shm
```

```bash
df | sed -n -r '2s/(..)%/\1/p'
#/dev/sda1            147386176  16617400 123161200  12 /
```

```bash
df | sed -n -r '2s/.*(..)%.*/\1/p'
#12
```

```bash
/sbin/ifconfig eth1
#eth1      Link encap:Ethernet  HWaddr 00:1B:11:15:E5:D8
#          inet addr:10.44.204.213  Bcast:10.44.207.255  Mask:255.255.252.0
#          inet6 addr: fe80::21b:11ff:fe15:e5d8/64 Scope:Link
#          UP BROADCAST RUNNING MULTICAST  MTU:1500  Metric:1
#          RX packets:5910584 errors:0 dropped:0 overruns:0 frame:0
#          TX packets:1141732 errors:0 dropped:0 overruns:0 carrier:0
#          collisions:0 txqueuelen:1000
#          RX bytes:463268584 (441.8 MiB)  TX bytes:224900586 (214.4 MiB)
#          Interrupt:217 Base address:0xef00
#
```

```bash
/sbin/ifconfig eth1 | sed -n -r '2s/.*:(.*)  (.*)  (.*)/\1/p'
#10.44.204.213
```

```bash
/sbin/ifconfig eth1 | sed -n -r '2s/^[^:]*:([^ ]*).*/\1/p'		#better way
#10.44.204.213
```

## for loop

```bash
vim for-0.bash
chmod +x for-0.bash

cat for-0.bash
##!/bin/bash
#
#for var in shantanu avinash deepak
#do
#        echo "Hello ${var}"
#done
```

```bash
./for-0.bash
#Hello shantanu
#Hello avinash
#Hello deepak
```

```bash
vim for-args.bash
chmod +x for-args.bash

cat for-args.bash
##!/bin/bash
#
#for var in $*
#do
#        echo "Hello ${var}"
#done

./for-args.bash

./for-args.bash shantanu avinash deepak
#Hello shantanu
#Hello avinash
#Hello deepak
```

```bash
who | cut -d " " -f1
#bashshell18
#bashshell3
#bashshell11
#bashshell1
#bashshell10
#bashshell4
#bashshell15
#bashshell25
#bashshell19
#bashshell19
#bashshell26
#bashshell27
#bashshell22
#bashshell21
#bashshell12
#bashshell12
#bashshell3
#bashshell20
#bashshell17
#bashshell8
#bashshell5
#bashshell9
#bashshell2
#bashshell24
#bashshell14
#bashshell28
```

```bash
who | cut -d " " -f1 | tr '\n' ' '
#bashshell18 bashshell3 bashshell11 bashshell1 bashshell10 bashshell4 bashshell15 bashshell25 bashshell19 bashshell19 bashshell26 bashshell27 bashshell22 bashshell21 bashshell12 bashshell12 bashshell3 bashshell20 bashshell17 bashshell8 bashshell5 bashshell9 bashshell2 bashshell24 bashshell14 bashshell28

vim for-who.bash

cat for-who.bash
##!/bin/bash
#
#for var in `who | cut -d ' ' -f1`
#do
#        echo "Hello ${var}"
#done

./for-who.bash
#Hello bashshell18
#Hello bashshell3
#Hello bashshell11
#Hello bashshell1
#Hello bashshell10
#Hello bashshell4
#Hello bashshell15
#Hello bashshell25
#Hello bashshell19
#Hello bashshell19
#Hello bashshell26
#Hello bashshell27
#Hello bashshell22
#Hello bashshell21
#Hello bashshell12
#Hello bashshell12
#Hello bashshell3
#Hello bashshell20
#Hello bashshell17
#Hello bashshell8
#Hello bashshell5
#Hello bashshell9
#Hello bashshell2
#Hello bashshell24
#Hello bashshell14
#Hello bashshell28
```

```bash
vim for-fileCount.bash

cat for-fileCount.bash
##!/bin/bash
#
#fc=0
#dc=0
#
#for var in *
#do
#        [[ -f ${var} ]] && let fc=fc+1
#        [[ -d ${var} ]] && let dc=dc+1
#        #echo ${var}
#done
#
#echo "$fc $dc"

./for-fileCount.bash
#28 0
```

```bash
vim for-fetchDataFromFile.bash

cat for-fetchDataFromFile.bash
##!/bin/bash
#
#cd /tmp/emails || exit 1
#
#for i in *@*
#do
#        (
#        echo $i
#        cat "$i"
#        ) | tr '\n' ' '
#        echo
#done

./for-fetchDataFromFile.bash
#ankita_srivasatva@persistent.com Ankita Srivastava
#awinashgautam@gmail.com awinash gautam
#awinashgautam@gmail.comcd
#hemalathav01@gmail.com Hemalatha V
#monali_khadke@persistent.com Monali Khadke
#nagarkar.pushkar@gmail.com Pushkar Nagarkar
#pratik_deshmukh@persistent.co.in
#pratik.jogwar@gmail.com Pratik Jogwar
#psbarbade@gmail.com
#rucha_hingne@persistent.co.in rucha
#cat: rucha_hingne@persistent.com: Is a directory
#rucha_hingne@persistent.com
#rupesh.pasalkar@gmail.com Rupesh Pasalkar
#sachin_vaze@persistent.co.in
#sachin_vpatil@persistent.co.in Sachin Patil
#sayali_gupte@persistent.co.in Sayali Gupte
#shahrohit2989@gmail.com Rohit Shah
#shahrohit2989@hotmail.com Rohit Shah
#sonali.hatapaki@persistent.co.in Sonali Hatapaki
#soniya_chavan@persistent.co.in Soniya
#tejasnjagtap@gmail.com tejas jagtap
```

```bash
vim while-1.bash

cat while-1.bash
##!/bin/bash
#
#while read a b c
#do
#        echo $a
#done < d.txt

./while-1.bash
#In
#A
#Where
#Through
#Down
```

```bash
vim while-1.bash

cat while-1.bash
##!/bin/bash
#
#IFS=','
#while read cname pop cap cpop junk
#do
#        echo "$cname    $cap"
#done < c.txt

./while-1.bash
#United Kingdom  London
#United States   Washington DC
#Venezuela       Caracas
#Vietnam Hanoi
#Yemen   San\'a
#Argentina       Buenos Aires
#Brazil  Brasilia
#Bahrain Manama
#Cameroon        Yaounde
#Djibouti        Djibouti
#Equatorial Guinea       Malabo
#Fiji    Suva
#France  Paris
#Greece  Athens
#Germany Berlin
#Honduras        Tegucigalpa
#China   Beijing
#Canada  Ottawa
#Hungary Budapest
#India   New Delhi
#Italy   Rome
#Ireland Dublin
#Japan   Tokio
```

```bash
vim while-2.bash

cat while-2.bash
##!/bin/bash
#
#who | cut -d ' ' -f1 | while read user
#do
#        echo "Hello $user"
#done

./while-2.bash
#Hello bashshell18
#Hello bashshell3
#Hello bashshell11
#Hello bashshell1
#Hello bashshell10
#Hello bashshell4
#Hello bashshell15
#Hello bashshell25
#Hello bashshell19
#Hello bashshell19
#Hello bashshell26
#Hello bashshell27
#Hello bashshell22
#Hello bashshell21
#Hello bashshell12
#Hello bashshell12
#Hello bashshell3
#Hello bashshell20
#Hello bashshell17
#Hello bashshell8
#Hello bashshell5
#Hello bashshell9
#Hello bashshell2
#Hello bashshell24
#Hello bashshell14
#Hello bashshell28
#Hello bashshell11
#Hello bashshell10
```

## functions

```bash
vim fun1.bash

chmod +x fun1.bash

cat fun1.bash
##!/bin/bash
#
#foo(){
#        echo "$1 -- $2 -- $3"
#        echo $#
#}
#
#foo hi bye ok globe

./fun1.bash
#hi -- bye -- ok
#4
```

```bash
vim fun1.bash

cat fun1.bash
##!/bin/bash
#
#foo(){
#        echo "$1 -- $2 -- $3"
#        echo $#
#        return 8        #return own custom value
#}
#
#foo hi bye ok globe
#echo "Return val: $?"

./fun1.bash
#hi -- bye -- ok
#4
#Return val: 8
```

```bash
vim fun2.bash

cat fun2.bash
##!/bin/bash
#
#a=10
#b=20
#
#foo(){
#        a=11            #a variable declared inside a function is global
#        local b         #use local keyword to make a var local
#        b=30
#}
#
#echo "Before calling function: a=${a}  b=${b}"
#foo
#echo "After calling function: a=${a}  b=${b}"
#

chmod +x fun2.bash

./fun2.bash
#Before calling function: a=10  b=20
#After calling function: a=11  b=20
```

```bash
vim dollar.bash
chmod +x dollar.bash

cat dollar.bash
##!/bin/bash
#
#for i in $*
#do
#        echo "hello $i"
#done

./dollar.bash foo bar car
#hello foo
#hello bar
#hello car
```

```bash
vim dollar.bash

cat dollar.bash
##!/bin/bash
#
#for i in "$*"
#do
#        echo "hello $i"
#done

./dollar.bash foo bar car
#hello foo bar car
```

```bash
vim dollar.bash

cat dollar.bash
##!/bin/bash
#
#for i in $@
#do
#        echo "hello $i"
#done

./dollar.bash foo bar car
#hello foo
#hello bar
#hello car
```

```bash
vim dollar.bash

cat dollar.bash
##!/bin/bash
#
#for i in "$@"
#do
#        echo "hello $i"
#done

./dollar.bash foo bar car
#hello foo
#hello bar
#hello car
```

#### Rule6: \$* and $@ are same except when they are in " "

```bash
# for "$@": it takes all the arguments as separate double quoted strings
# for "$*": it takes all the arguments as a single double quoted string where all the arguments are separated by space
```


## process management


```bash
sleep 10000 &
#[1] 11838

jobs
#[1]+  Running                 sleep 10000 &

sleep 10002 &
#[2] 11841

jobs
#[1]-  Running                 sleep 10000 &
#[2]+  Running                 sleep 10002 &
```

```bash
sleep 10005 &
#[3] 11878

jobs
#[1]   Running                 sleep 10000 &
#[2]-  Running                 sleep 10002 &
#[3]+  Running                 sleep 10005 &

sleep 10007 &
#[4] 11910

sleep 10009 &
#[5] 11911

jobs
#[1]   Running                 sleep 10000 &
#[2]   Running                 sleep 10002 &
#[3]   Running                 sleep 10005 &
#[4]-  Running                 sleep 10007 &
#[5]+  Running                 sleep 10009 &

jobs
#[1]   Running                 sleep 10000 &
#[2]   Running                 sleep 10002 &
#[3]   Running                 sleep 10005 &
#[4]-  Running                 sleep 10007 &
#[5]+  Running                 sleep 10009 &

jobs
#[1]   Running                 sleep 10000 &
#[2]   Running                 sleep 10002 &
#[3]   Running                 sleep 10005 &
#[4]-  Running                 sleep 10007 &
#[5]+  Running                 sleep 10009 &

#get the process to foreground
fg %2        
#sleep 10002
#
#[2]+  Stopped                 sleep 10002

#after pressing Ctrl+z, the process comes stops

#then resume the process in background by bg %2

bg %2
#[2]+ sleep 10002 &

jobs
#[1]   Running                 sleep 10000 &
#[2]   Running                 sleep 10002 &
#[3]   Running                 sleep 10005 &
#[4]-  Running                 sleep 10007 &
#[5]+  Running                 sleep 10009 &
```

## sort

```bash
cat c.txt
#United Kingdom,57533000,London,6756000,Europe,1066,1945,-,English
#United States,252177000,Washington DC,606900,North America,1776,1945,-,English
#Venezuela,19733000,Caracas,1290087,Latin America,1811,1945,-,Spanish
#Vietnam,67589000,Hanoi,1088862,Orient,1945,1977,-,Vietnamese
#Yemen,1184300,San\'a,427150,Asia,1918,1957,Islam,Arabic
#Argentina,36955182,Buenos Aires,2033445,Latin America,1853,1945,-,Spanish
#Brazil,172860370,Brasilia,286037,Latin America,1822,1945,-,Portuguese
#Bahrain,634137,Manama,34137,Persian Gulf,1973,1977,Islamic,Arabic
#Cameroon,15421937,Yaounde,421937,Africa,1960,1974,-,Franch
#Djibouti,451442,Djibouti,1442,Africa,1977,1980,-,Franch
#Equatorial Guinea,474214,Malabo,74214,Africa,1991,1995,-,Franch
#Fiji,832494,Suva,32494,Oceania,1970,1975,-,English
#France,59329691,Paris,329691,Europe,486,1945,-,Franch
#Greece,10601527,Athens,601527,Europe,1829,1945,-,Greek
#Germany,82797408,Berlin,1797408,Europe,1871,1960,-,German
#Honduras,6249598,Tegucigalpa,1249598,Latin America,1821,1945,-,Spanish
#China,1261832482,Beijing,61832482,Asia,-221,1945,-,Chinese
#Canada,31281092,Ottawa,1281092,North America,1867,1945,-,English
#Hungary,10138844,Budapest,1138844,Europe,1001,1945,-,Hungerian
#India,1014003817,New Delhi,14003817,Asia,1947,1950,-,Indian
#Italy,57634327,Rome,3634327,Europe,1861,1950,-,Italian
#Ireland,3797257,Dublin,797257,Europe,1921,1945,-,English
#Japan,126549976,Tokio,16549976,Asia,-660,1955,-,Japanese
```

```bash
sort -t, -k3 c.txt
#Greece,10601527,Athens,601527,Europe,1829,1945,-,Greek
#China,1261832482,Beijing,61832482,Asia,-221,1945,-,Chinese
#Germany,82797408,Berlin,1797408,Europe,1871,1960,-,German
#Brazil,172860370,Brasilia,286037,Latin America,1822,1945,-,Portuguese
#Hungary,10138844,Budapest,1138844,Europe,1001,1945,-,Hungerian
#Argentina,36955182,Buenos Aires,2033445,Latin America,1853,1945,-,Spanish
#Venezuela,19733000,Caracas,1290087,Latin America,1811,1945,-,Spanish
#Djibouti,451442,Djibouti,1442,Africa,1977,1980,-,Franch
#Ireland,3797257,Dublin,797257,Europe,1921,1945,-,English
#Vietnam,67589000,Hanoi,1088862,Orient,1945,1977,-,Vietnamese
#United Kingdom,57533000,London,6756000,Europe,1066,1945,-,English
#Equatorial Guinea,474214,Malabo,74214,Africa,1991,1995,-,Franch
#Bahrain,634137,Manama,34137,Persian Gulf,1973,1977,Islamic,Arabic
#India,1014003817,New Delhi,14003817,Asia,1947,1950,-,Indian
#Canada,31281092,Ottawa,1281092,North America,1867,1945,-,English
#France,59329691,Paris,329691,Europe,486,1945,-,Franch
#Italy,57634327,Rome,3634327,Europe,1861,1950,-,Italian
#Yemen,1184300,San\'a,427150,Asia,1918,1957,Islam,Arabic
#Fiji,832494,Suva,32494,Oceania,1970,1975,-,English
#Honduras,6249598,Tegucigalpa,1249598,Latin America,1821,1945,-,Spanish
#Japan,126549976,Tokio,16549976,Asia,-660,1955,-,Japanese
#United States,252177000,Washington DC,606900,North America,1776,1945,-,English
#Cameroon,15421937,Yaounde,421937,Africa,1960,1974,-,Franch
```

```bash
#sort by 9th col and stop there. next col to consider for sort is 3rd col
sort -t, -k9,9 -k3 c.txt     
#Bahrain,634137,Manama,34137,Persian Gulf,1973,1977,Islamic,Arabic
#Yemen,1184300,San\'a,427150,Asia,1918,1957,Islam,Arabic
#China,1261832482,Beijing,61832482,Asia,-221,1945,-,Chinese
#Ireland,3797257,Dublin,797257,Europe,1921,1945,-,English
#United Kingdom,57533000,London,6756000,Europe,1066,1945,-,English
#Canada,31281092,Ottawa,1281092,North America,1867,1945,-,English
#Fiji,832494,Suva,32494,Oceania,1970,1975,-,English
#United States,252177000,Washington DC,606900,North America,1776,1945,-,English
#Djibouti,451442,Djibouti,1442,Africa,1977,1980,-,Franch
#Equatorial Guinea,474214,Malabo,74214,Africa,1991,1995,-,Franch
#France,59329691,Paris,329691,Europe,486,1945,-,Franch
#Cameroon,15421937,Yaounde,421937,Africa,1960,1974,-,Franch
#Germany,82797408,Berlin,1797408,Europe,1871,1960,-,German
#Greece,10601527,Athens,601527,Europe,1829,1945,-,Greek
#Hungary,10138844,Budapest,1138844,Europe,1001,1945,-,Hungerian
#India,1014003817,New Delhi,14003817,Asia,1947,1950,-,Indian
#Italy,57634327,Rome,3634327,Europe,1861,1950,-,Italian
#Japan,126549976,Tokio,16549976,Asia,-660,1955,-,Japanese
#Brazil,172860370,Brasilia,286037,Latin America,1822,1945,-,Portuguese
#Argentina,36955182,Buenos Aires,2033445,Latin America,1853,1945,-,Spanish
#Venezuela,19733000,Caracas,1290087,Latin America,1811,1945,-,Spanish
#Honduras,6249598,Tegucigalpa,1249598,Latin America,1821,1945,-,Spanish
#Vietnam,67589000,Hanoi,1088862,Orient,1945,1977,-,Vietnamese
```

```bash
#sort by 9th col and 3rd col 2nd char
sort -t, -k9,9 -k3.2 c.txt | grep English    
#United States,252177000,Washington DC,606900,North America,1776,1945,-,English
#United Kingdom,57533000,London,6756000,Europe,1066,1945,-,English
#Canada,31281092,Ottawa,1281092,North America,1867,1945,-,English
#Ireland,3797257,Dublin,797257,Europe,1921,1945,-,English
#Fiji,832494,Suva,32494,Oceania,1970,1975,-,English
```

```bash
#uniq
cut -d, -f9 c.txt | sort -u          
#Arabic
#Chinese
#English
#Franch
#German
#Greek
#Hungerian
#Indian
#Italian
#Japanese
#Portuguese
#Spanish
#Vietnamese
```

```bash
#uniq command
cut -d, -f9 c.txt | sort | uniq              
#Arabic
#Chinese
#English
#Franch
#German
#Greek
#Hungerian
#Indian
#Italian
#Japanese
#Portuguese
#Spanish
#Vietnamese
```

```bash
#uniq command with -c which is not available with sort -u
cut -d, -f9 c.txt | sort | uniq -c           
#      2 Arabic
#      1 Chinese
#      5 English
#      4 Franch
#      1 German
#      1 Greek
#      1 Hungerian
#      1 Indian
#      1 Italian
#      1 Japanese
#      1 Portuguese
#      3 Spanish
#      1 Vietnamese
```

```bash
#awk: alphred peter_wineburger brian_kernignen
```

```bash
cat -A e.txt
#SR^I8649^I275^IAsia$
#Canada^I3852^I25^INorth America$
#China^I3705^I1032^IAsia$
#USA^I3615^I237^INorth America$
#Brazil^I3286^I134^ISouth America$
#India^I1267^I746^IAsia$
#Mexico^I762^I78^INorth America$
#France^I211^I55^IEurope$
#Japan^I144^I120^IAsia$
#Germany^I96^I61^IEurope$
#England^I94^I56^IEurope$
```

## awk 'Pattern {action}' filename

```bash
#awk opens file to read line by line
#pattern is searched on each line and the action is performed
#either Pattern or action can be omitted but not both
#if pattern is omitted, action is performed on each line
#if action is omitted, a default action of print is performed on the line in which pattern is found
```

```bash
awk '{ print "hello", $0 }' e.txt
#hello SR        8649    275     Asia
#hello Canada    3852    25      North America
#hello China     3705    1032    Asia
#hello USA       3615    237     North America
#hello Brazil    3286    134     South America
#hello India     1267    746     Asia
#hello Mexico    762     78      North America
#hello France    211     55      Europe
#hello Japan     144     120     Asia
#hello Germany   96      61      Europe
#hello England   94      56      Europe
```

```bash
awk '{ print "hello", $1 }' e.txt
#hello SR
#hello Canada
#hello China
#hello USA
#hello Brazil
#hello India
#hello Mexico
#hello France
#hello Japan
#hello Germany
#hello England
```

```bash
awk '{ print "country", $1, "continent", $4 }' e.txt
#country SR continent Asia
#country Canada continent North
#country China continent Asia
#country USA continent North
#country Brazil continent South
#country India continent Asia
#country Mexico continent North
#country France continent Europe
#country Japan continent Asia
#country Germany continent Europe
#country England continent Europe
```

```bash
awk -F"\t" '{ print "country", $1, "continent", $4 }' e.txt
#country SR continent Asia
#country Canada continent North America
#country China continent Asia
#country USA continent North America
#country Brazil continent South America
#country India continent Asia
#country Mexico continent North America
#country France continent Europe
#country Japan continent Asia
#country Germany continent Europe
#country England continent Europe
```

```bash
#prints number of fields
awk -F"\t" '{ print NF }' e.txt      
#4
#4
#4
#4
#4
#4
#4
#4
#4
#4
#4
```

```bash
#last column: NF is field/column number
awk -F"\t" '{ print $NF }' e.txt     
#Asia
#North America
#Asia
#North America
#South America
#Asia
#North America
#Europe
#Asia
#Europe
#Europe
```

```bash
#second last column
awk -F"\t" '{ print $(NF-1) }' e.txt         
#275
#25
#1032
#237
#134
#746
#78
#55
#120
#61
#56
```

```bash
#second last column * 10
awk -F"\t" '{ print $(NF-1)*10 }' e.txt      
#2750
#250
#10320
#2370
#1340
#7460
#780
#550
#1200
#610
#560
```

```bash
#NR is row number
awk '{ print NR, $0 }' e.txt         
#1 SR    8649    275     Asia
#2 Canada        3852    25      North America
#3 China 3705    1032    Asia
#4 USA   3615    237     North America
#5 Brazil        3286    134     South America
#6 India 1267    746     Asia
#7 Mexico        762     78      North America
#8 France        211     55      Europe
#9 Japan 144     120     Asia
#10 Germany      96      61      Europe
#11 England      94      56      Europe
```

```bash
awk '/Asia/' e.txt
#SR      8649    275     Asia
#China   3705    1032    Asia
#India   1267    746     Asia
#Japan   144     120     Asia
```

```bash
#~ : contains
awk '$NF ~ /America/' e.txt          
#Canada  3852    25      North America
#USA     3615    237     North America
#Brazil  3286    134     South America
#Mexico  762     78      North America
```

```bash
#~ : contains
awk -F"\t" '$NF ~ /^N.*ca$/' e.txt           
#Canada  3852    25      North America
#USA     3615    237     North America
#Mexico  762     78      North America
```

```bash
awk -F"\t" '$NF == "Asia"' e.txt
#SR      8649    275     Asia
#China   3705    1032    Asia
#India   1267    746     Asia
#Japan   144     120     Asia
```

```bash
awk -F"\t" '$2 > 3000' e.txt
#SR      8649    275     Asia
#Canada  3852    25      North America
#China   3705    1032    Asia
#USA     3615    237     North America
#Brazil  3286    134     South America
```

```bash
awk 'NR>9 && NR<16 {print NR,$0}' c.txt
#10 Djibouti,451442,Djibouti,1442,Africa,1977,1980,-,Franch
#11 Equatorial Guinea,474214,Malabo,74214,Africa,1991,1995,-,Franch
#12 Fiji,832494,Suva,32494,Oceania,1970,1975,-,English
#13 France,59329691,Paris,329691,Europe,486,1945,-,Franch
#14 Greece,10601527,Athens,601527,Europe,1829,1945,-,Greek
#15 Germany,82797408,Berlin,1797408,Europe,1871,1960,-,German
```

```bash
cal
#     June 2015
#Su Mo Tu We Th Fr Sa
#    1  2  3  4  5  6
# 7  8  9 10 11 12 13
#14 15 16 17 18 19 20
#21 22 23 24 25 26 27
#28 29 30
#

#task: show 15 from this cal

cal | awk -F" " 'NR==5 {print $2}'
#15
```

```bash
df
#Filesystem           1K-blocks      Used Available Use% Mounted on
#/dev/sda1            147386176  16636604 123141996  12% /
#tmpfs                  1032324         0   1032324   0% /dev/shm

#task: show 12

df | awk 'NR==2 { print int($5) }'
#12

df | awk 'NR==2 { printf "%d\n", $5 }'
#12

df | awk 'NR==2 { sub("%", "", $5); print $5 }'
#12
```

```bash
/sbin/ifconfig eth1
#eth1      Link encap:Ethernet  HWaddr 00:1B:11:15:E5:D8
#          inet addr:10.44.204.213  Bcast:10.44.207.255  Mask:255.255.252.0
#          inet6 addr: fe80::21b:11ff:fe15:e5d8/64 Scope:Link
#          UP BROADCAST RUNNING MULTICAST  MTU:1500  Metric:1
#          RX packets:11092960 errors:0 dropped:0 overruns:0 frame:0
#          TX packets:5115509 errors:0 dropped:0 overruns:0 carrier:0
#          collisions:0 txqueuelen:1000
#          RX bytes:927938064 (884.9 MiB)  TX bytes:490855054 (468.1 MiB)
#          Interrupt:217 Base address:0xef00
#

#task: sho the ip address

/sbin/ifconfig eth1 | awk 'NR==2 { print $3 }'
#Bcast:10.44.207.255

/sbin/ifconfig eth1 | awk 'NR==2 { print $2 }'
#addr:10.44.204.213

/sbin/ifconfig eth1 | awk 'NR==2 { sub("addr:", "", $2); print $2 }'
#10.44.204.213
```

```bash
cal
#     June 2015
#Su Mo Tu We Th Fr Sa
#    1  2  3  4  5  6
# 7  8  9 10 11 12 13
#14 15 16 17 18 19 20
#21 22 23 24 25 26 27
#28 29 30
#

#task: show 15 and 22 both on same line

cal | awk 'NR==5 || NR==6 {print $2}'
#15
#22

cal | awk 'NR==5 || NR==6 { sub("\n", " ", $2); print $2}'
#15
#22

cal | awk 'NR==5 || NR==6 { printf("%d ", $2)}'
#15 22

cal | awk 'NR==5 { printf("%d ", $2)} NR==6 {print $2}'
#15 22

cal | awk 'NR==5 { var=$2 } NR==6 {print var, $2}'
#15 22
```

```bash
cat e.txt
#SR      8649    275     Asia
#Canada  3852    25      North America
#China   3705    1032    Asia
#USA     3615    237     North America
#Brazil  3286    134     South America
#India   1267    746     Asia
#Mexico  762     78      North America
#France  211     55      Europe
#Japan   144     120     Asia
#Germany 96      61      Europe
#England 94      56      Europe

#task: print the name of countries whose continent is Europe

awk -F"\t" '$NF=="Europe" {print $1}' e.txt
#France
#Germany
#England
```

```bash
#BEGIN block is executed before 1st line is read
awk -F"\t" 'BEGIN { print "List of countries: " }
#> $NF=="Europe" { print $1 }' e.txt
#List of countries:
#France
#Germany
#England
```

```bash
#END block is executed after file is read
awk -F"\t" 'BEGIN { print "List of countries: " }
#$NF=="Europe" { print $1 }
#> END { print "Thank you" }' e.txt
#List of countries:
#France
#Germany
#England
#Thank you
```

```bash
awk -F"\t" 'BEGIN { c=0 }
#$NF=="Europe" { c++ }
#END { print "European countries " c}' e.txt
#European countries 3


awk -F"\t" '$NF=="Europe" { c++ }
#END { print "European countries " c}' e.txt
#European countries 3


awk -F"\t" '$NF=="Europe" { c=c+$3 }
#END { print "Total population in European countries " c}' e.txt
#Total population in European countries 172


awk -F"\t" '$NF=="Europe" { c=c " " $1 }
#END { print "All European countries: " c}' e.txt
#All European countries:  France Germany England
```


#### associative arrays


```bash
awk 'BEGIN { X["a"]="bar"; X["b"]="foo"; X["c"]="car"; print X["a"] }'
#bar

awk 'BEGIN { X["a"]="bar"; X["b"]="foo"; X["c"]="car"; i="a"; print X[i] }'
#bar

awk 'BEGIN { X["a"]="bar"; X["b"]="foo"; X["c"]="car";
#for(i in X){ print i, X[i] } }'
#a bar
#b foo
#c car
```

```bash
#associative arrays to count the frequency
cat e.txt
#SR      8649    275     Asia
#Canada  3852    25      North America
#China   3705    1032    Asia
#USA     3615    237     North America
#Brazil  3286    134     South America
#India   1267    746     Asia
#Mexico  762     78      North America
#France  211     55      Europe
#Japan   144     120     Asia
#Germany 96      61      Europe
#England 94      56      Europe
```

```bash
awk -F"\t" '{ X[$NF]++ }
#END{ for(i in X){ print i, X[i]}}' e.txt
#Asia 4
#Europe 3
#South America 1
#North America 3
```

```bash
awk -F"\t" '{ X[$NF]=X[$NF] " " $1}
#END{ for(i in X){ print i, " ---> ", X[i]}}' e.txt
#Asia  --->   SR China India Japan
#Europe  --->   France Germany England
#South America  --->   Brazil
#North America  --->   Canada USA Mexico

awk -F"\t" '{ X[$NF]=X[$NF] " " $1}
#END{ for(i in X){ print i, " ---> ", X[i]}}' e.txt
#Asia  --->   SR China India Japan
#Europe  --->   France Germany England
#South America  --->   Brazil
#North America  --->   Canada USA Mexico
```

```bash
vim f1

cat f1
#sneha
#rags
#abhya
#snehal
#madhura

vim f2

cat f2
#rewati
#rags
#bhagi
#madhura
#shruti

sort f1 > f3
sort f2 > f4

cat f3
#abhya
#madhura
#rags
#sneha
#snehal

cat f4
#bhagi
#madhura
#rags
#rewati
#shruti

comm f3 f4
#abhya
#        bhagi
#                madhura
#                rags
#        rewati
#        shruti
#sneha
#snehal

#1st col: lines uniq to f3
#2nd col: lines uniq to f4
#3rd col: lines uniq to f3 and f4

comm -12 f3 f4       #shows only 3rd col
#madhura
#rags
```

```bash
#rename all *.html files to *.txt file
vim fRename.bash

cat fRename.bash
##!/bin/bash
#
#for file in *.html
#do
#        mv "$file" "${file%html}txt"
#done
```

```bash
./fRename.bash
```

```bash
vim case1.bash
chmod +x case1.bash

cat case1.bash
##!/bin/bash
#
#echo 'Enter name: '
#read name
#
#case "$name" in
#        david) echo 'David Korn'
#        ;;
#        bill) echo 'Bill Joy'
#        ;;
#        chet|brian) echo 'Authors of bash'
#        ;;
#        *) echo 'Who is this?'
#esac
```

```bash
./case1.bash
#Enter name:
#bill
#Bill Joy
```

## Books

```bash
#Learning the Bash Shell: Cameroon Newham
#Classic Shell Scripting: Arnold Robins
#Shell Scripting for Automation: Arnold Robins
#Sed and Awk: Arnold Robins
#Learning the vi editor: Arnold Robins
#Unix Power Tools: Best book to become expert
```
