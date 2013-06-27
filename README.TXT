###BIOJAZZ


####INSTALLATION


Perl and certain CPAN modules must be installed on your system.

#####Perl Modules installation

I have found the following to work when installing
perl modules (e.g. Bit::Vector) as a user.
```
cd ~
    rm -rf .cpan
    export PERL5LIB=~/myperl  # (and add to .bashrc)
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
The project is firstly developed by Julien Ollivier, then modified and currently maintained by Song Feng (@LifeWorks).

### Support or Contact
If have any problems, please contact @LifeWorks or email lifeworks.song@gmail.com
