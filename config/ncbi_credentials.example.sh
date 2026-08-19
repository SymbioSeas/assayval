# AssayVal NCBI credentials — OPTIONAL, for faster assembly downloads
#
# Set your NCBI API key here ONCE and the download helper
# (scripts/download/download_assemblies.sh) applies it on every run.
#
# Setup:
#   cp config/ncbi_credentials.example.sh config/ncbi_credentials.sh
#   # then edit config/ncbi_credentials.sh and paste your key below
#
# config/ncbi_credentials.sh is gitignored, so your key is never committed.
#
# Get a free key: https://www.ncbi.nlm.nih.gov/account/ (Account Settings →
# "API Key Management"). A key raises your NCBI rate limit from 3 to 10
# requests/sec, which meaningfully speeds up large multi-thousand-genome pulls.
#
# Precedence: a -k flag on the command line overrides an NCBI_API_KEY environment
# variable, which overrides the value set in this file.

export NCBI_API_KEY="paste-your-key-here"

# Optional. The `datasets` CLI does not require an e-mail, but you can set one
# here for your own records / other NCBI tooling.
# export NCBI_EMAIL="you@example.org"
