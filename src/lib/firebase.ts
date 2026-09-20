import { initializeApp, getApps, getApp } from 'firebase/app';
import { getAuth } from 'firebase/auth';
import { getFirestore } from 'firebase/firestore';

// Helper to safely read env vars in both Vite client and Node environments
const getEnv = (key: string): string => {
  try {
    if (typeof process !== 'undefined' && process.env && process.env[key]) {
      return process.env[key] as string;
    }
    // @ts-ignore
    if (typeof import.meta !== 'undefined' && import.meta.env && import.meta.env[key]) {
      // @ts-ignore
      return import.meta.env[key] as string;
    }
  } catch {
    // ignore
  }
  return '';
};

// Safe configuration using environment variables or resilient fallback defaults
const firebaseConfig = {
  apiKey: getEnv('VITE_FIREBASE_API_KEY') || 'AIzaSyDemoInteractionDegradationApiKey',
  authDomain: getEnv('VITE_FIREBASE_AUTH_DOMAIN') || 'interaction-app.firebaseapp.com',
  projectId: getEnv('VITE_FIREBASE_PROJECT_ID') || 'interaction-app',
  storageBucket: getEnv('VITE_FIREBASE_STORAGE_BUCKET') || 'interaction-app.appspot.com',
  messagingSenderId: getEnv('VITE_FIREBASE_MESSAGING_SENDER_ID') || '100000000000',
  appId: getEnv('VITE_FIREBASE_APP_ID') || '1:100000000000:web:abcdef1234567890',
  firestoreDatabaseId: getEnv('VITE_FIREBASE_FIRESTORE_DATABASE_ID') || '(default)'
};

const app = getApps().length === 0 ? initializeApp(firebaseConfig) : getApp();
export const db = getFirestore(app, firebaseConfig.firestoreDatabaseId);
export const auth = getAuth(app);
export default app;
