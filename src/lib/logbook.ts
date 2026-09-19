import { db } from './firebase';
import { 
  collection, 
  addDoc, 
  getDocs, 
  query, 
  orderBy, 
  limit, 
  doc, 
  setDoc, 
  getDoc, 
  deleteDoc, 
  writeBatch 
} from 'firebase/firestore';

export interface LogbookEntry {
  id?: string;
  primaryCompound: string;
  secondaryCompounds: string;
  predictedImpurities: string;
  timestamp: string;
  createdAt?: number;
}

/**
 * Escapes values according to RFC 4180 CSV specifications.
 */
export function escapeCsvCell(value: string | undefined | null): string {
  if (value === undefined || value === null) return '""';
  const str = String(value);
  return `"${str.replace(/"/g, '""')}"`;
}

/**
 * Extracts date portion (YYYY-MM-DD) from a timestamp string.
 */
export function formatDateOnly(value?: string): string {
  if (!value) return '';
  const str = String(value).trim();
  if (str.includes('T')) return str.split('T')[0];
  if (str.includes(' ')) return str.split(' ')[0];
  return str;
}

/**
 * Builds standard CSV content from an array of logbook entries.
 * Columns: primary compound, secondary compounds, predicted impurities, timestamp (date only)
 */
export function buildCsvString(entries: LogbookEntry[]): string {
  const header = ['"primary_compound"', '"secondary_compounds"', '"predicted_impurities"', '"timestamp"'].join(',');
  const rows = entries.map(e => [
    escapeCsvCell(e.primaryCompound),
    escapeCsvCell(e.secondaryCompounds || 'None'),
    escapeCsvCell(e.predictedImpurities || 'None'),
    escapeCsvCell(formatDateOnly(e.timestamp))
  ].join(','));

  return [header, ...rows].join('\n');
}

/**
 * Records a chemical interaction query to Firestore, maintaining only the 100 most recent queries,
 * and updates the consolidated CSV logbook document on the database.
 */
export async function logQueryToDatabase(
  primarySmiles: string,
  secondarySmilesList: string[],
  impurities: Array<{ name?: string; iupacName?: string; smiles?: string; probability?: number }>
): Promise<void> {
  try {
    const timestamp = new Date().toISOString().split('T')[0];
    const createdAt = Date.now();
    const secondaryStr = secondarySmilesList.filter(s => s && s.trim()).join('; ') || 'None';

    // Format predicted impurities (name + SMILES)
    const impuritiesStr = impurities && impurities.length > 0
      ? impurities.map(imp => {
          const name = imp.iupacName || imp.name || 'Unknown Product';
          const smiles = imp.smiles ? ` [${imp.smiles}]` : '';
          const prob = imp.probability != null ? ` (${(imp.probability * 100).toFixed(1)}%)` : '';
          return `${name}${smiles}${prob}`;
        }).join('; ')
      : 'None detected';

    const newEntry: LogbookEntry = {
      primaryCompound: primarySmiles.trim(),
      secondaryCompounds: secondaryStr,
      predictedImpurities: impuritiesStr,
      timestamp,
      createdAt
    };

    // 1. Add new query log
    await addDoc(collection(db, 'query_logs'), newEntry);

    // 2. Fetch all query logs ordered by createdAt desc to maintain the 100-limit window
    const logsQuery = query(collection(db, 'query_logs'), orderBy('createdAt', 'desc'));
    const snapshot = await getDocs(logsQuery);
    const docs = snapshot.docs;

    // Prune entries older than 100
    if (docs.length > 100) {
      const docsToDelete = docs.slice(100);
      const batch = writeBatch(db);
      docsToDelete.forEach(d => batch.delete(d.ref));
      await batch.commit();
    }

    // Keep top 100 entries for the CSV logbook
    const recentDocs = docs.slice(0, 100);
    const recentEntries: LogbookEntry[] = recentDocs.map(d => ({
      id: d.id,
      ...(d.data() as LogbookEntry)
    }));

    // 3. Update the CSV logbook document on the database
    const csvContent = buildCsvString(recentEntries);
    await setDoc(doc(db, 'logbook', 'csv_logbook'), {
      csvContent,
      entryCount: recentEntries.length,
      updatedAt: timestamp
    });

  } catch (error) {
    console.error('Failed to log query to database:', error);
  }
}

/**
 * Retrieves the current list of up to 100 recent query logs from the database.
 */
export async function getRecentQueryLogs(maxResults: number = 100): Promise<LogbookEntry[]> {
  try {
    const q = query(
      collection(db, 'query_logs'), 
      orderBy('createdAt', 'desc'), 
      limit(maxResults)
    );
    const snap = await getDocs(q);
    return snap.docs.map(d => ({ id: d.id, ...(d.data() as LogbookEntry) }));
  } catch (err) {
    console.error('Error fetching query logs:', err);
    return [];
  }
}

/**
 * Retrieves the maintained CSV logbook from Firestore.
 */
export async function getCsvLogbookFromDatabase(): Promise<string> {
  try {
    const logbookDoc = await getDoc(doc(db, 'logbook', 'csv_logbook'));
    if (logbookDoc.exists() && logbookDoc.data().csvContent) {
      return logbookDoc.data().csvContent;
    }
    // Fallback: build on the fly from query_logs if csv_logbook doc not created yet
    const entries = await getRecentQueryLogs(100);
    return buildCsvString(entries);
  } catch (err) {
    console.error('Error fetching CSV logbook from database:', err);
    return '"primary_compound","secondary_compounds","predicted_impurities","timestamp"\n';
  }
}

/**
 * Purges all collections from the database.
 */
export async function wipeAllDatabaseData(): Promise<{ success: boolean; error?: string }> {
  try {
    const collectionsToPurge = ['query_logs', 'logbook', 'compounds', 'predictions'];
    for (const coll of collectionsToPurge) {
      const snap = await getDocs(collection(db, coll));
      if (!snap.empty) {
        const batch = writeBatch(db);
        snap.docs.forEach(d => batch.delete(d.ref));
        await batch.commit();
      }
    }
    return { success: true };
  } catch (err: any) {
    console.error('Error wiping database:', err);
    return { success: false, error: err.message || 'Wipe operation failed' };
  }
}
