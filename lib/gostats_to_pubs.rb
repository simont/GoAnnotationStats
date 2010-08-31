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
unique_pmids = Hash.new

Bio::GO::GeneAssociation.parser(go_data) do |entry|
  # puts "Found Ref: #{entry.db_reference.to_s}"
  pmids = entry.db_reference.to_s.scan(/PMID\:(\d+)/)
  if !pmids.empty?
    pmids.each do |pmid|
      unique_pmids.has_key?(pmid) ? unique_pmids[pmid] += 1 : unique_pmids[pmid] = 1
    end
  end
end

puts "Found #{unique_pmids.size} unique PMIDs in the GAF list"
# unique_pmids.each do |k,v|
#   puts "#{k} has #{v.to_s} annotation"
# end
