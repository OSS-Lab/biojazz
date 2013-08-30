###BIOJAZZ


####INSTALLATION


Perl and certain CPAN modules must be installed on your system.

#####Perl Modules installation

I have found the following to work when installing
perl modules (e.g. Bit::Vector) as a user.
```
cd ~

perl -MCPAN -e shell

o conf make_arg -I/users/songfeng/CPAN

o conf make_install_arg -I/users/songfeng/CPAN

o conf makepl_arg "install_base=/users/songfeng/CPAN LIB=/users/songfeng/CPAN/lib INSTALLPRIVLIB=/users/songfeng/CPAN/lib INSTALLARCHLIB=/users/songfeng/CPAN/lib/arch INSTALLSITEARCH=/users/songfeng/CPAN/lib/arch INSTALLSITELIB=/users/songfeng/CPAN/lib INSTALLSCRIPT=/users/songfeng/CPAN/bin INSTALLBIN=/users/songfeng/CPAN/bin INSTALLSITEBIN=/users/songfeng/CPAN/bin INSTALLMAN1DIR=/users/songfeng/CPAN/man/man1 INSTALLSITEMAN1DIR=/users/songfeng/CPAN/man/man1 INSTALLMAN3DIR=/users/songfeng/CPAN/man/man3 INSTALLSITEMAN3DIR=/users/songfeng/CPAN/man/man3"

    export PERL5LIB=/users/songfeng/CPAN  # (and add to .bashrc)
perl -MCPAN -e shell   # (or just cpan if in path)
```
answer the questions, and when asked, say LIB=~/myperl PREFIX=~/myperl
```
    cpan> install Bit::Vector
```

####USAGE

    Use the -h option to get usage info:

    E.g. biojazz.pl -h




### Authors and Contributors
The project is firstly developed by Julien Ollivier, then modified and currently maintained by @LifeWorks.

### Support or Contact
If have any problems, please contact @LifeWorks.
