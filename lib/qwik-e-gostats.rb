#########
# Simon Twigger, RGD, 2009-2010
#
# Script to parse GO Consortium Gene Annotation files and count basic annotation stats
#
# Uses Bio::Ruby to provide the GO Annotation file parsing
#
#########

require 'rubygems'
require 'bio'
require 'digest/md5'

EVIDENCE_CODES = ['EXP','IDA','IPI','IMP','IGI','IEP','ISS','ISO','ISA','ISM','IGC','RCA','TAS','NAS','IC','ND','IEA']
MANUAL_CODES = ['EXP','IDA','IPI','IMP','IGI','IEP','IGC','TAS','NAS','IC','ND']

file = ARGV[0] || '../spec/files/gene_association_short.rgd'

go_data = File.open(file).read
go = Bio::GO::GeneAssociation.parser(go_data)

source = Hash.new
gene_products = Hash.new
unique_annotations = Hash.new

Bio::GO::GeneAssociation.parser(go_data) do |entry|
  str = Array.new
  str   << entry.db << entry.db_object_id << entry.db_object_symbol << entry.goid << entry.with.to_s << entry.qualifier  << entry.evidence << entry.with.to_s << entry.db_object_name << entry.db_object_type << entry.db_reference.to_s << entry.taxon << entry.assigned_by

  digest = Digest::MD5.hexdigest(str.to_s)
  unique_annotations.has_key?(digest) ? unique_annotations[digest]['count'] += 1 : unique_annotations[digest] = {'count'  => 1, 'str'  => str}
  gene_products.has_key?(entry.entry_id) ? gene_products[entry.entry_id] += 1 : gene_products[entry.entry_id] = 1
  if !source.has_key?(entry.assigned_by)
    source[entry.assigned_by] = { 'products' => {entry.entry_id => 1}, 'annotations' => {
      'EXP' => 0,
      'IDA' => 0,
      'IPI' => 0,
      'IMP' => 0,
      'IGI' => 0,
      'IEP' => 0,
      'ISS' => 0,
      'ISO' => 0,
      'ISA' => 0,
      'ISM' => 0,
      'IGC' => 0,
      'RCA' => 0,
      'TAS' => 0,
      'NAS' => 0,
      'IC' => 0,
      'ND' => 0,
      'IEA' => 0,
      },
      'manual_annotations' => {}}
    end

    # source[entry.assigned_by]['annotations']['total'] += 1
    source[entry.assigned_by]['products'].has_key?(entry.entry_id) ? source[entry.assigned_by]['products'][entry.entry_id] += 1 : source[entry.assigned_by]['products'][entry.entry_id] = 1
    source[entry.assigned_by]['annotations'][entry.evidence] += 1
    # if entry.evidence != "IEA" && entry.evidence != "ISO" && (entry.evidence != "ISS" && entry.assigned_by == "RGD") # ignore ISS for RGD's manual count
    if MANUAL_CODES.include?(entry.evidence)
      source[entry.assigned_by]['manual_annotations'].has_key?(entry.entry_id) ? source[entry.assigned_by]['manual_annotations'][entry.entry_id] += 1 : source[entry.assigned_by]['manual_annotations'][entry.entry_id] = 1
    end

  end




  puts "Unique gene products:\t#{gene_products.size}"
  puts "Unique sources:\t#{source.size}"
  # puts "Unique annotation rows: #{unique_annotations.size}"

  # unique_annotations.sort_by{ |annot| annot[1][:count]}.each do |k|
  total_rows = 0
  total_dups = 0
  unique_annotations.each do |k,v|

    if unique_annotations[k]['count'] > 1
      # puts unique_annotations[k]['count'].to_s << " copies of: " << unique_annotations[k]['str'].join(' ')
      total_dups += unique_annotations[k]['count']-1
    end
    total_rows += unique_annotations[k]['count']
  end

  puts "Total number of annotation rows:\t#{total_rows}"
  puts "Total number of non-unique annotation rows:\t#{total_dups}"
  puts "Total number of unique annotation rows:\t#{total_rows - total_dups}"

  sorted_sources = source.keys.sort

  # source.each_pair do |group,data|
  sorted_sources.each do |group|
    total = source[group]['annotations'].each_value.inject {|total, ec| total + ec}
    source[group]['total_annotations'] = total
  end
  print "EC\t"
  puts sorted_sources.join("\t")

  EVIDENCE_CODES.each do |ec|
    ec_str = "#{ec}"
    # put an asterisk after manual codes
    MANUAL_CODES.include?(ec) ? ec_str << "*\t" : ec_str << "\t"

    print ec_str
    ec_total = 0
    # source.keys.each do |group|
    sorted_sources.each do |group|
      ec_total += source[group]['annotations'][ec]
      print "#{source[group]['annotations'][ec]}\t"
    end
    print "#{ec_total}\n"

  end

  print "Total Annotations:\t"
  # source.each_pair do |group,data|
  sorted_sources.each do |group|
    print "#{source[group]['total_annotations']}\t"
  end
  print "\n"

  print "Unique Genes Annotated:\t"
  # source.each_pair do |group,data|
  sorted_sources.each do |group|
    print "#{source[group]['products'].size}\t"
  end
  print "\n"

  print "Total Manual Genes Annotated:\t"
  sorted_sources.each do |group|
    print "#{source[group]['manual_annotations'].size}\t"
  end
  print "\n"
  
  print "Total Manual Annotations:\t"
  sorted_sources.each do |group|
    print "#{source[group]['manual_annotations'].each_value.inject {|sum,tot| sum + tot}}\t"
  end
  print "\n"
