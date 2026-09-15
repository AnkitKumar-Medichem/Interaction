import React, { useState, useEffect } from 'react';
import { 
  FileSpreadsheet, 
  Download, 
  Copy, 
  Check, 
  RefreshCw, 
  Trash2, 
  Search, 
  Database,
  Calendar,
  AlertCircle
} from 'lucide-react';
import { 
  LogbookEntry, 
  getRecentQueryLogs, 
  getCsvLogbookFromDatabase, 
  wipeAllDatabaseData,
  buildCsvString 
} from '../lib/logbook';

interface CsvLogbookProps {
  onNewReactionClick?: () => void;
  onSelectSmiles?: (smiles: string) => void;
}

export const CsvLogbook: React.FC<CsvLogbookProps> = ({ onNewReactionClick, onSelectSmiles }) => {
  const [entries, setEntries] = useState<LogbookEntry[]>([]);
  const [csvRaw, setCsvRaw] = useState<string>('');
  const [loading, setLoading] = useState<boolean>(true);
  const [searchTerm, setSearchTerm] = useState<string>('');
  const [copied, setCopied] = useState<boolean>(false);
  const [activeTab, setActiveTab] = useState<'table' | 'raw'>('table');
  const [showWipeConfirm, setShowWipeConfirm] = useState<boolean>(false);
  const [wiping, setWiping] = useState<boolean>(false);
  const [statusMessage, setStatusMessage] = useState<{ text: string; type: 'success' | 'error' } | null>(null);

  const fetchLogbook = async () => {
    setLoading(true);
    try {
      const logs = await getRecentQueryLogs(100);
      setEntries(logs);
      const raw = await getCsvLogbookFromDatabase();
      setCsvRaw(raw || buildCsvString(logs));
    } catch (err) {
      console.error('Failed to load logbook:', err);
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    fetchLogbook();
  }, []);

  const handleCopyCsv = async () => {
    try {
      const textToCopy = csvRaw || buildCsvString(entries);
      await navigator.clipboard.writeText(textToCopy);
      setCopied(true);
      setTimeout(() => setCopied(false), 2000);
    } catch (err) {
      console.error('Failed to copy to clipboard', err);
    }
  };

  const handleDownloadCsv = () => {
    const textToDownload = csvRaw || buildCsvString(entries);
    const blob = new Blob([textToDownload], { type: 'text/csv;charset=utf-8;' });
    const url = URL.createObjectURL(blob);
    const link = document.createElement('a');
    const timestamp = new Date().toISOString().replace(/[:.]/g, '-').slice(0, 19);
    link.setAttribute('href', url);
    link.setAttribute('download', `interaction_logbook_100_${timestamp}.csv`);
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    URL.revokeObjectURL(url);
  };

  const handleWipeDatabase = async () => {
    setWiping(true);
    try {
      const res = await wipeAllDatabaseData();
      if (res.success) {
        setEntries([]);
        setCsvRaw('"primary_compound","secondary_compounds","predicted_impurities","timestamp"\n');
        setStatusMessage({ text: 'Database successfully cleared. All query records and logbook wiped.', type: 'success' });
        setShowWipeConfirm(false);
        setTimeout(() => setStatusMessage(null), 4000);
      } else {
        setStatusMessage({ text: res.error || 'Failed to wipe database.', type: 'error' });
      }
    } catch (err: any) {
      setStatusMessage({ text: err.message || 'Error executing database wipe.', type: 'error' });
    } finally {
      setWiping(false);
    }
  };

  const filteredEntries = entries.filter(e => {
    if (!searchTerm) return true;
    const term = searchTerm.toLowerCase();
    return (
      e.primaryCompound.toLowerCase().includes(term) ||
      (e.secondaryCompounds && e.secondaryCompounds.toLowerCase().includes(term)) ||
      (e.predictedImpurities && e.predictedImpurities.toLowerCase().includes(term)) ||
      e.timestamp.toLowerCase().includes(term)
    );
  });

  return (
    <div className="w-full space-y-6">
      {/* Top Banner & Controls */}
      <div className="bg-white border border-[#E2E8F0] rounded-2xl p-6 sm:p-8 shadow-[0_1px_3px_rgba(15,23,42,0.03)]">
        <div className="flex flex-col md:flex-row md:items-center justify-between gap-4 pb-6 border-b border-[#F1F5F9]">
          <div>
            <div className="flex items-center gap-2 mb-1">
              <FileSpreadsheet className="w-6 h-6 text-[#4F46E5]" />
              <h2 className="font-serif text-2xl font-bold text-[#0F172A]">Database CSV Logbook</h2>
              <span className="text-xs font-semibold px-2.5 py-0.5 rounded-full bg-[#EEF2FF] text-[#4F46E5] border border-[#E0E7FF]">
                {entries.length} / 100 Recent Queries
              </span>
            </div>
            <p className="text-sm text-[#64748B]">
              Persistently maintained on the database. Contains the most recent 100 SMILES interaction queries with columns: primary compound, secondary compounds, predicted impurities, and timestamp.
            </p>
          </div>

          <div className="flex flex-wrap items-center gap-2">
            <button
              onClick={fetchLogbook}
              disabled={loading}
              className="inline-flex items-center gap-1.5 px-3 py-2 text-xs font-semibold text-[#334155] bg-[#F8FAFC] hover:bg-[#F1F5F9] border border-[#E2E8F0] rounded-lg transition-colors"
              title="Refresh database logbook"
            >
              <RefreshCw className={`w-3.5 h-3.5 ${loading ? 'animate-spin' : ''}`} />
              Refresh
            </button>
            <button
              onClick={handleCopyCsv}
              disabled={entries.length === 0}
              className="inline-flex items-center gap-1.5 px-3 py-2 text-xs font-semibold text-[#334155] bg-[#F8FAFC] hover:bg-[#F1F5F9] border border-[#E2E8F0] rounded-lg transition-colors disabled:opacity-50"
            >
              {copied ? <Check className="w-3.5 h-3.5 text-green-600" /> : <Copy className="w-3.5 h-3.5" />}
              {copied ? 'Copied CSV' : 'Copy CSV'}
            </button>
            <button
              onClick={handleDownloadCsv}
              disabled={entries.length === 0}
              className="inline-flex items-center gap-1.5 px-4 py-2 text-xs font-semibold text-white bg-[#4F46E5] hover:bg-[#4338CA] rounded-lg transition-colors shadow-xs disabled:opacity-50"
            >
              <Download className="w-3.5 h-3.5" />
              Download CSV Logbook
            </button>
            <button
              onClick={() => setShowWipeConfirm(true)}
              className="inline-flex items-center gap-1.5 px-3 py-2 text-xs font-semibold text-red-600 hover:text-red-700 bg-red-50 hover:bg-red-100 border border-red-200 rounded-lg transition-colors"
              title="Delete everything from database"
            >
              <Trash2 className="w-3.5 h-3.5" />
              Clear Database
            </button>
          </div>
        </div>

        {statusMessage && (
          <div className={`mt-4 p-3 rounded-lg text-xs font-medium flex items-center gap-2 ${
            statusMessage.type === 'success' ? 'bg-green-50 text-green-800 border border-green-200' : 'bg-red-50 text-red-800 border border-red-200'
          }`}>
            <AlertCircle className="w-4 h-4 shrink-0" />
            {statusMessage.text}
          </div>
        )}

        {/* Confirmation Modal / Banner for Wipe */}
        {showWipeConfirm && (
          <div className="mt-4 p-4 bg-red-50/80 border border-red-200 rounded-xl">
            <div className="text-sm font-bold text-red-900 mb-1">Confirm Complete Database Wipe</div>
            <p className="text-xs text-red-700 mb-3">
              This will irreversibly delete all stored query logs and the CSV logbook document from the Firestore database.
            </p>
            <div className="flex items-center gap-2">
              <button
                onClick={handleWipeDatabase}
                disabled={wiping}
                className="px-3 py-1.5 bg-red-600 hover:bg-red-700 text-white text-xs font-semibold rounded-lg transition-colors disabled:opacity-50"
              >
                {wiping ? 'Wiping...' : 'Yes, Delete Everything'}
              </button>
              <button
                onClick={() => setShowWipeConfirm(false)}
                disabled={wiping}
                className="px-3 py-1.5 bg-white border border-[#CBD5E1] text-[#334155] hover:bg-[#F8FAFC] text-xs font-semibold rounded-lg transition-colors"
              >
                Cancel
              </button>
            </div>
          </div>
        )}

        {/* Search and View Mode Tabs */}
        <div className="mt-6 flex flex-col sm:flex-row items-center justify-between gap-3">
          <div className="relative w-full sm:w-80">
            <Search className="w-4 h-4 text-[#94A3B8] absolute left-3 top-1/2 -translate-y-1/2" />
            <input
              type="text"
              value={searchTerm}
              onChange={(e) => setSearchTerm(e.target.value)}
              placeholder="Search by SMILES or predicted impurity..."
              className="w-full pl-9 pr-3 py-2 text-xs bg-[#F8FAFC] border border-[#E2E8F0] rounded-lg text-[#0F172A] placeholder:text-[#94A3B8] focus:outline-none focus:border-[#4F46E5] focus:ring-1 focus:ring-[#4F46E5]"
            />
          </div>

          <div className="flex items-center bg-[#F1F5F9] p-1 rounded-lg self-end sm:self-auto">
            <button
              onClick={() => setActiveTab('table')}
              className={`px-3 py-1 text-xs font-semibold rounded-md transition-colors ${
                activeTab === 'table' ? 'bg-white text-[#0F172A] shadow-xs' : 'text-[#64748B] hover:text-[#0F172A]'
              }`}
            >
              Table View ({filteredEntries.length})
            </button>
            <button
              onClick={() => setActiveTab('raw')}
              className={`px-3 py-1 text-xs font-semibold rounded-md transition-colors ${
                activeTab === 'raw' ? 'bg-white text-[#0F172A] shadow-xs' : 'text-[#64748B] hover:text-[#0F172A]'
              }`}
            >
              Raw CSV Output
            </button>
          </div>
        </div>

        {/* Content Views */}
        <div className="mt-4">
          {loading ? (
            <div className="py-16 text-center">
              <div className="inline-block w-8 h-8 border-3 border-[#EEF2FF] border-t-[#4F46E5] rounded-full animate-spin mb-3"></div>
              <p className="text-xs text-[#64748B]">Retrieving logbook from Firestore...</p>
            </div>
          ) : entries.length === 0 ? (
            <div className="py-16 text-center bg-[#F8FAFC] border border-dashed border-[#CBD5E1] rounded-xl">
              <Database className="w-10 h-10 text-[#94A3B8] mx-auto mb-2" />
              <div className="text-sm font-semibold text-[#0F172A]">Database Logbook is Empty</div>
              <p className="text-xs text-[#64748B] max-w-sm mx-auto mt-1 mb-4">
                No chemical reaction queries logged yet. Run a prediction using SMILES inputs to automatically record it in the database logbook.
              </p>
              {onNewReactionClick && (
                <button
                  onClick={onNewReactionClick}
                  className="px-4 py-2 bg-[#4F46E5] hover:bg-[#4338CA] text-white text-xs font-semibold rounded-lg shadow-xs transition-colors"
                >
                  Start Reaction Setup
                </button>
              )}
            </div>
          ) : activeTab === 'table' ? (
            <div className="overflow-x-auto border border-[#E2E8F0] rounded-xl">
              <table className="w-full text-left text-xs text-[#334155]">
                <thead className="bg-[#F8FAFC] border-b border-[#E2E8F0] text-[11px] font-bold text-[#475569] uppercase tracking-wider">
                  <tr>
                    <th className="py-3 px-4 w-12 text-center">#</th>
                    <th className="py-3 px-4 min-w-[200px]">Primary Compound (SMILES)</th>
                    <th className="py-3 px-4 min-w-[180px]">Secondary Compounds</th>
                    <th className="py-3 px-4 min-w-[280px]">Predicted Impurities</th>
                    <th className="py-3 px-4 min-w-[160px]">Timestamp</th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-[#F1F5F9]">
                  {filteredEntries.map((entry, idx) => (
                    <tr key={entry.id || idx} className="hover:bg-[#F8FAFC] transition-colors">
                      <td className="py-3 px-4 text-center font-mono text-[#94A3B8]">
                        {idx + 1}
                      </td>
                      <td className="py-3 px-4">
                        <div className="font-mono text-[11px] bg-slate-50 border border-slate-200 px-2 py-1 rounded text-slate-800 break-all select-all font-medium">
                          {entry.primaryCompound}
                        </div>
                      </td>
                      <td className="py-3 px-4">
                        {entry.secondaryCompounds && entry.secondaryCompounds !== 'None' ? (
                          <div className="font-mono text-[11px] bg-indigo-50/50 border border-indigo-100 px-2 py-1 rounded text-indigo-900 break-all select-all">
                            {entry.secondaryCompounds}
                          </div>
                        ) : (
                          <span className="text-slate-400 italic">None</span>
                        )}
                      </td>
                      <td className="py-3 px-4 text-[#0F172A]">
                        <div className="max-h-24 overflow-y-auto leading-relaxed text-[11px] whitespace-pre-wrap">
                          {entry.predictedImpurities || 'None detected'}
                        </div>
                      </td>
                      <td className="py-3 px-4 text-[#64748B] whitespace-nowrap">
                        <div className="flex items-center gap-1.5 font-mono text-[11px]">
                          <Calendar className="w-3.5 h-3.5 text-[#94A3B8]" />
                          {entry.timestamp ? new Date(entry.timestamp).toLocaleString('en-US', {
                            year: 'numeric',
                            month: 'short',
                            day: 'numeric',
                            hour: '2-digit',
                            minute: '2-digit',
                            second: '2-digit',
                            hour12: false
                          }) : 'N/A'}
                        </div>
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          ) : (
            <div className="relative">
              <textarea
                readOnly
                value={csvRaw || buildCsvString(entries)}
                rows={14}
                className="w-full p-4 font-mono text-xs bg-[#0F172A] text-emerald-400 rounded-xl border border-slate-700 leading-relaxed focus:outline-none"
              />
              <button
                onClick={handleCopyCsv}
                className="absolute top-3 right-3 px-3 py-1.5 bg-slate-800 hover:bg-slate-700 text-slate-200 rounded-lg text-xs font-medium transition-colors border border-slate-600 flex items-center gap-1.5"
              >
                {copied ? <Check className="w-3 h-3 text-emerald-400" /> : <Copy className="w-3 h-3" />}
                {copied ? 'Copied' : 'Copy'}
              </button>
            </div>
          )}
        </div>
      </div>
    </div>
  );
};
