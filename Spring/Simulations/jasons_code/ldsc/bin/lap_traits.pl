BEGIN{
    #push @INC,'/broad/software/free/Linux/redhat_6_x86_64/pkgs/perl_5.10.1/lib/5.10.1/';
    #push @INC,'/broad/software/free/Linux/redhat_6_x86_64/pkgs/perl_5.10.1/lib/site_perl/5.10.1/';
    push @INC,'/home/unix/tgreen/perl5/lib/site_perl/5.10.1/';
    #push @INC,'/home/unix/tgreen/dig-jenkins/lib/';
}


use strict;
use Getopt::Long;
use Pod::Usage;
use Config::Simple;
use DBI;
use Log::Log4perl;
#use KeyArgs;

my $help;
my $man;
my $dbi_configuration_file;
my $debug;
my $trace;
my $autocommit = 1;
my $log_file_name='dump_property.log';
my @param = @ARGV;
my $database_table;
my $no_phenotype;
my $ancestry;
my $no_ancestry;
my $parent;
my @avoid_datasets = ();
GetOptions(
    'trace'                         => \$trace,
    'debug'                         => \$debug,
    'help'                          => \$help,
    'man'                           => \$man,
    'database-table=s'              => \$database_table,
    'no-phenotype'                  => \$no_phenotype,
    'ancestry=s'                      => \$ancestry,
    'no-ancestry=s'                      => \$no_ancestry,
    'no-dataset=s'                      => \@avoid_datasets,
    'parent=s'                      => \$parent,
    'dbi-configuration-file=s'      => \$dbi_configuration_file,
    ) or pod2usage(-verbose         => 1) && exit;

pod2usage({-noperldoc => 1, -verbose => 1,-exitval=>1}) if defined $help;
pod2usage({-noperldoc => 1, -verbose => 2,-exitval=>1}) if defined $man;
pod2usage({-message=>"\n\n ERROR: --dbi-configuration-file not configured\n\n",-noperldoc => 1, -verbose => 2}) if (! -e $dbi_configuration_file );

my $log4perlconf = qq(
    log4perl.category.Foo.Bar        = INFO, Screen , Log
    log4perl.category.debugger       = DEBUG, Screen
    log4perl.category.tracer         = TRACE, Screen
    log4perl.appender.Screen         = Log::Log4perl::Appender::Screen
    log4perl.appender.Screen.stderr  = 1
    log4perl.appender.Screen.layout= PatternLayout
    log4perl.appender.Screen.layout.ConversionPattern=  %d %p %m%n
    log4perl.appender.Log            = Log::Log4perl::Appender::File
    log4perl.appender.Log.filename   = $log_file_name
    log4perl.appender.Log.mode       = append
    log4perl.appender.Log.layout= PatternLayout
    log4perl.appender.Log.layout.ConversionPattern=  %d %p %m%n
  );

Log::Log4perl::init( \$log4perlconf );
my $logger;
if($debug){
     $logger = Log::Log4perl->get_logger('debugger');
}elsif($trace){
     $logger = Log::Log4perl->get_logger('tracer');
}else{
     $logger = Log::Log4perl->get_logger('Foo.Bar');
}

$logger->info("Starting $0 @param");

my $cfg;
my $db_string;
my $host_url;
my $db_password;
my $dbh;
my $user;

$cfg = new Config::Simple($dbi_configuration_file);
$user=$cfg->param('user');
$db_string=$cfg->param('database');
$host_url=$cfg->param('host_url');
$db_password=$cfg->param('db_pswd');
$dbh =  DBI->connect("DBI:mysql:database=${db_string};host=${host_url}",$user , $db_password, 
		     {'RaiseError' => 1, 'AutoCommit' => $autocommit});

my $avoid = "";
if (@avoid_datasets)
{
		$avoid = " AND " . join(" AND ", map {/(.*),(.*)/ ? "(n1.PH != '$1' || n1.DATASET != '$2')" : "n1.DATASET != '$_'"} @avoid_datasets);
}

sub get_dataset_phenotype_table($$){
		
    my @return;
    #my $sql = "SELECT  META_MDV.DATASET,PH,TBL FROM
    #           META_MDV,META_DATASET_PH
    #           WHERE META_MDV.DATASET=META_DATASET_PH.DATASET";

		my $sql = "select n2.DATASET, n1.PH, n1.TBL, CASES, SUBJECTS, n2.ANCESTRY from META_DATASET_PH n1 INNER JOIN META_DATASET n2 INNER JOIN (select PH, MAX(SUBJECTS) m from META_DATASET_PH n1, META_DATASET m2 WHERE n1.DATASET = m2.DATASET and m2.EXPTYPE = 'GWAS' $avoid GROUP BY PH) b where b.PH = n1.PH and n1.DATASET = n2.DATASET and b.m = n2.SUBJECTS $avoid";

		if ($ancestry)
		{
				$sql = "select n2.DATASET, n1.PH, n1.TBL, CASES, SUBJECTS, n2.ANCESTRY from META_DATASET_PH n1 INNER JOIN META_DATASET n2 INNER JOIN (select PH, MAX(SUBJECTS) m from META_DATASET_PH n1, META_DATASET m2 WHERE n1.DATASET = m2.DATASET and m2.EXPTYPE = 'GWAS' $avoid GROUP BY PH) b where b.PH = n1.PH and n1.DATASET = n2.DATASET and b.m = n2.SUBJECTS AND n2.ANCESTRY = '$ancestry' $avoid";
		}
		elsif ($no_ancestry) 
		{
				$sql = "select n2.DATASET, n1.PH, n1.TBL, CASES, SUBJECTS, n2.ANCESTRY from META_DATASET_PH n1 INNER JOIN META_DATASET n2 INNER JOIN (select PH, MAX(SUBJECTS) m from META_DATASET_PH n1, META_DATASET m2 WHERE n1.DATASET = m2.DATASET and m2.EXPTYPE = 'GWAS' $avoid GROUP BY PH) b where b.PH = n1.PH and n1.DATASET = n2.DATASET and b.m = n2.SUBJECTS AND n2.ANCESTRY != '$no_ancestry' $avoid";
		} 

		my $dbh = shift;
		my $logger = shift;
		my $exp_filter = shift;

    if($exp_filter)
		{
				$sql .= "\nAND EXP like '$exp_filter'"
    }

		print STDERR "$sql\n";
    $logger->debug($sql);
    my $sth = $dbh->prepare($sql);
    $sth->execute();
    while(my @this = $sth->fetchrow()){
				push @return,{DATASET=>$this[0],
											PH=>$this[1],
											TBL=>$this[2],
											CASES=>$this[3],
											SUBJECTS=>$this[4],
											ANCESTRY=>$this[5]}
    }
    $sth->finish();
    return \@return;
}

my $tbl_phenotype = get_dataset_phenotype_table($dbh,
																								$logger);
my %seen;
foreach my $row (@$tbl_phenotype){
    next if $seen{$row->{TBL}};
    print "$row->{TBL}\t$row->{PH}\t$row->{DATASET}\t$row->{CASES}\t$row->{SUBJECTS}\t$row->{ANCESTRY}\n";
    $seen{$row->{TBL}}=1;
}
$dbh->disconnect();

=head1 NAME

 lap_dataset.pl

=head1 SYNOPSIS

 creates dataset

=head1 DESCRIPTION

 

=head1 OPTIONS

 debug                       set logger to DEBUG
 trace                       set logger to TRACE
 help                        help
 man                         manual

 dbi-configuration-file


=head1 EXAMPLE

 perl lap_dataset.pl --dbi-config /humgen/diabetes/portal/continuous_integration/servers/sandbox/resources/aws_dev.dbi
 
=cut
