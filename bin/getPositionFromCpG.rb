require 'selenium-webdriver'
require 'csv'

cpg_csv = CSV.read("./data/coef.csv")

caps = Selenium::WebDriver::Remote::Capabilities.chrome("chromeOptions" => {args: ["--headless", "--disable-gpu"]})
driver = Selenium::WebDriver.for :chrome, desired_capabilities: caps

url = "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19" 

driver.navigate.to url

cpg_position = cpg_csv[2..-1].map do |cpg|
  form = driver.find_element(:name, 'hgt.positionInput')
  form.send_keys(cpg[1])
  driver.find_element(:name, "hgt.jump").click

  sleep(1)

  driver.find_element(:id, "positionDisplay").text
end

CSV.open( "./data/tsd_cpgs.csv", "w" ){ |csv|
  csv << ["cpg id", "coef", "chr", "position"]
  csv << [cpg_csv[1][1], cpg_csv[1][2], "", ""]
  cpg_csv[2..-1].each_with_index do |cpg, i|
    position = cpg_position[i].split(":")
    chr = position[0]
    position = position[1].split("-")[0].gsub(/(\d{0,3}),(\d{3})/, '\1\2').to_i
    csv << [cpg[1], cpg[2], chr, position]
  end
}


