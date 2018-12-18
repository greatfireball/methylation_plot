#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

use Log::Log4perl;

# Initialize Logger
my $log_conf = q(
   log4perl.rootLogger              = DEBUG, SCREEN
   log4perl.appender.SCREEN         = Log::Log4perl::Appender::Screen
   log4perl.appender.SCREEN.stderr  = 0
   log4perl.appender.SCREEN.layout  = Log::Log4perl::Layout::PatternLayout
   log4perl.appender.SCREEN.layout.ConversionPattern = [%d] %p %m %n
);
Log::Log4perl::init(\$log_conf);
my $log = Log::Log4perl->get_logger();
